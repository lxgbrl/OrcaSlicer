#include "ImplicitExpr.hpp"

#include <cctype>
#include <cmath>
#include <memory>
#include <stdexcept>

namespace Slic3r {

namespace {
struct ParseError : std::runtime_error { using std::runtime_error::runtime_error; };
constexpr double kPi = 3.14159265358979323846;
} // namespace

// Recursive-descent parser. Compiles each node directly into a std::function
// closure (Fn). Symbol resolution uses the owning ImplicitExpr so that calls
// bind to user functions / constants known at this point in the program.
struct ImplicitExpr::Parser
{
    const std::string& s;
    std::size_t        i = 0;
    ImplicitExpr&      owner;

    Parser(const std::string& str, ImplicitExpr& o) : s(str), owner(o) {}

    void skip() { while (i < s.size() && std::isspace((unsigned char)s[i])) ++i; }
    bool eof()  { skip(); return i >= s.size(); }
    char peek() { skip(); return i < s.size() ? s[i] : '\0'; }
    bool accept(char c) { if (peek() == c) { ++i; return true; } return false; }
    bool accept2(char a, char b) { skip(); if (i + 1 < s.size() && s[i] == a && s[i+1] == b) { i += 2; return true; } return false; }

    [[noreturn]] void fail(const std::string& m) { throw ParseError(m + " at pos " + std::to_string(i)); }

    Fn parse_top()
    {
        Fn f = parse_or();
        if (!eof()) fail("unexpected trailing input");
        return f;
    }

    // || (lowest)
    Fn parse_or()
    {
        Fn a = parse_and();
        for (;;) {
            if (accept2('|', '|')) { Fn b = parse_and(); a = [a, b](const Env& e){ return (a(e) != 0.0 || b(e) != 0.0) ? 1.0 : 0.0; }; }
            else return a;
        }
    }
    Fn parse_and()
    {
        Fn a = parse_cmp();
        for (;;) {
            if (accept2('&', '&')) { Fn b = parse_cmp(); a = [a, b](const Env& e){ return (a(e) != 0.0 && b(e) != 0.0) ? 1.0 : 0.0; }; }
            else return a;
        }
    }
    Fn parse_cmp()
    {
        Fn a = parse_add();
        for (;;) {
            skip();
            if (accept2('<', '=')) { Fn b = parse_add(); a = [a,b](const Env& e){ return a(e) <= b(e) ? 1.0 : 0.0; }; }
            else if (accept2('>', '=')) { Fn b = parse_add(); a = [a,b](const Env& e){ return a(e) >= b(e) ? 1.0 : 0.0; }; }
            else if (accept2('=', '=')) { Fn b = parse_add(); a = [a,b](const Env& e){ return a(e) == b(e) ? 1.0 : 0.0; }; }
            else if (accept2('!', '=')) { Fn b = parse_add(); a = [a,b](const Env& e){ return a(e) != b(e) ? 1.0 : 0.0; }; }
            else if (peek() == '<')     { ++i; Fn b = parse_add(); a = [a,b](const Env& e){ return a(e) <  b(e) ? 1.0 : 0.0; }; }
            else if (peek() == '>')     { ++i; Fn b = parse_add(); a = [a,b](const Env& e){ return a(e) >  b(e) ? 1.0 : 0.0; }; }
            else return a;
        }
    }
    Fn parse_add()
    {
        Fn a = parse_mul();
        for (;;) {
            char c = peek();
            if (c == '+') { ++i; Fn b = parse_mul(); a = [a,b](const Env& e){ return a(e) + b(e); }; }
            else if (c == '-') { ++i; Fn b = parse_mul(); a = [a,b](const Env& e){ return a(e) - b(e); }; }
            else return a;
        }
    }
    Fn parse_mul()
    {
        Fn a = parse_unary();
        for (;;) {
            char c = peek();
            if (c == '*') { ++i; Fn b = parse_unary(); a = [a,b](const Env& e){ return a(e) * b(e); }; }
            else if (c == '/') { ++i; Fn b = parse_unary(); a = [a,b](const Env& e){ return a(e) / b(e); }; }
            else return a;
        }
    }
    Fn parse_unary()
    {
        char c = peek();
        if (c == '-') { ++i; Fn a = parse_unary(); return [a](const Env& e){ return -a(e); }; }
        if (c == '+') { ++i; return parse_unary(); }
        return parse_pow();
    }
    Fn parse_pow()
    {
        Fn a = parse_atom();
        skip();
        if (peek() == '^') { ++i; Fn b = parse_unary(); return [a,b](const Env& e){ return std::pow(a(e), b(e)); }; }
        return a;
    }

    Fn parse_number()
    {
        skip();
        std::size_t start = i;
        while (i < s.size() && (std::isdigit((unsigned char)s[i]) || s[i]=='.' )) ++i;
        if (i < s.size() && (s[i]=='e' || s[i]=='E')) { ++i; if (i<s.size() && (s[i]=='+'||s[i]=='-')) ++i; while (i<s.size() && std::isdigit((unsigned char)s[i])) ++i; }
        double v = std::stod(s.substr(start, i - start));
        return [v](const Env&){ return v; };
    }

    std::string parse_ident()
    {
        skip();
        std::size_t start = i;
        while (i < s.size() && (std::isalnum((unsigned char)s[i]) || s[i]=='_')) ++i;
        return s.substr(start, i - start);
    }

    std::vector<Fn> parse_args()
    {
        std::vector<Fn> args;
        if (!accept('(')) fail("expected '('");
        if (accept(')')) return args;
        for (;;) {
            args.push_back(parse_or());
            if (accept(',')) continue;
            if (accept(')')) break;
            fail("expected ',' or ')'");
        }
        return args;
    }

    static Fn builtin(const std::string& name, std::vector<Fn>& a)
    {
        auto need = [&](std::size_t n){ if (a.size() != n) throw ParseError("builtin '" + name + "' arity"); };
        if (name == "if") { need(3); Fn c=a[0],x=a[1],y=a[2]; return [c,x,y](const Env& e){ return c(e)!=0.0 ? x(e) : y(e); }; }
        if (name == "min") { if (a.size()<2) throw ParseError("min arity"); auto v=a; return [v](const Env& e){ double m=v[0](e); for (std::size_t k=1;k<v.size();++k) m=std::min(m,v[k](e)); return m; }; }
        if (name == "max") { if (a.size()<2) throw ParseError("max arity"); auto v=a; return [v](const Env& e){ double m=v[0](e); for (std::size_t k=1;k<v.size();++k) m=std::max(m,v[k](e)); return m; }; }
        if (name == "atan2") { need(2); Fn x=a[0],y=a[1]; return [x,y](const Env& e){ return std::atan2(x(e),y(e)); }; }
        if (name == "pow")   { need(2); Fn x=a[0],y=a[1]; return [x,y](const Env& e){ return std::pow(x(e),y(e)); }; }
        if (name == "mod" || name == "fmod") { need(2); Fn x=a[0],y=a[1]; return [x,y](const Env& e){ return std::fmod(x(e),y(e)); }; }
        if (a.size() == 1) {
            Fn x = a[0];
            if (name=="sin")  return [x](const Env& e){ return std::sin(x(e)); };
            if (name=="cos")  return [x](const Env& e){ return std::cos(x(e)); };
            if (name=="tan")  return [x](const Env& e){ return std::tan(x(e)); };
            if (name=="asin") return [x](const Env& e){ return std::asin(x(e)); };
            if (name=="acos") return [x](const Env& e){ return std::acos(x(e)); };
            if (name=="atan") return [x](const Env& e){ return std::atan(x(e)); };
            if (name=="sinh") return [x](const Env& e){ return std::sinh(x(e)); };
            if (name=="cosh") return [x](const Env& e){ return std::cosh(x(e)); };
            if (name=="tanh") return [x](const Env& e){ return std::tanh(x(e)); };
            if (name=="sqrt") return [x](const Env& e){ return std::sqrt(x(e)); };
            if (name=="cbrt") return [x](const Env& e){ return std::cbrt(x(e)); };
            if (name=="abs"||name=="fabs") return [x](const Env& e){ return std::fabs(x(e)); };
            if (name=="exp")  return [x](const Env& e){ return std::exp(x(e)); };
            if (name=="ln"||name=="log") return [x](const Env& e){ return std::log(x(e)); };
            if (name=="log10")return [x](const Env& e){ return std::log10(x(e)); };
            if (name=="sign") return [x](const Env& e){ double v=x(e); return (v>0)-(v<0)+0.0; };
            if (name=="floor")return [x](const Env& e){ return std::floor(x(e)); };
            if (name=="ceil") return [x](const Env& e){ return std::ceil(x(e)); };
            if (name=="round")return [x](const Env& e){ return std::round(x(e)); };
        }
        return nullptr; // not a builtin
    }

    // Build a call to a user function (def already compiled). Rebinds x,y,z,t.
    static Fn user_call(Fn callee, std::vector<Fn> args)
    {
        return [callee, args](const Env& e) {
            Env ne;
            ne.x = args.size() > 0 ? args[0](e) : e.x;
            ne.y = args.size() > 1 ? args[1](e) : e.y;
            ne.z = args.size() > 2 ? args[2](e) : e.z;
            ne.t = args.size() > 3 ? args[3](e) : e.t;
            return callee(ne);
        };
    }

    Fn parse_atom()
    {
        skip();
        char c = peek();
        if (c == '(') { ++i; Fn f = parse_or(); if (!accept(')')) fail("expected ')'"); return f; }
        if (std::isdigit((unsigned char)c) || c == '.') return parse_number();
        if (std::isalpha((unsigned char)c) || c == '_') {
            std::string name = parse_ident();
            skip();
            bool has_call = (peek() == '(');
            // bare variables / constants
            if (!has_call) {
                if (name=="x") return [](const Env& e){ return e.x; };
                if (name=="y") return [](const Env& e){ return e.y; };
                if (name=="z") return [](const Env& e){ return e.z; };
                if (name=="t") return [](const Env& e){ return e.t; };
                if (name=="pi"||name=="PI") { double v=kPi; return [v](const Env&){ return v; }; }
                auto cit = owner.m_consts.find(name);
                if (cit != owner.m_consts.end()) { double v = cit->second; return [v](const Env&){ return v; }; }
                auto fit = owner.m_funcIndex.find(name);
                if (fit != owner.m_funcIndex.end()) return owner.m_defs[fit->second]; // call with current env
                fail("unknown identifier '" + name + "'");
            }
            std::vector<Fn> args = parse_args();
            if (Fn b = builtin(name, args)) return b;
            auto fit = owner.m_funcIndex.find(name);
            if (fit != owner.m_funcIndex.end()) return user_call(owner.m_defs[fit->second], std::move(args));
            fail("unknown function '" + name + "'");
        }
        fail("unexpected character");
    }
};

bool ImplicitExpr::compile_formula(const std::string&              formula,
                                   const std::vector<std::string>& consts,
                                   std::string&                    error)
{
    return compile(consts, {}, formula, error);
}

bool ImplicitExpr::compile(const std::vector<std::string>& consts,
                           const std::vector<std::string>& funct,
                           const std::string&               fxyz,
                           std::string&                     error)
{
    m_root = nullptr;
    m_defs.clear();
    m_funcIndex.clear();
    m_consts.clear();
    try {
        // Constants: "name = expr" (expr may reference earlier constants).
        for (const std::string& line : consts) {
            auto eq = line.find('=');
            if (eq == std::string::npos) continue;
            std::string name = line.substr(0, eq);
            // trim
            name.erase(0, name.find_first_not_of(" \t"));
            name.erase(name.find_last_not_of(" \t") + 1);
            std::string rhs = line.substr(eq + 1);
            Parser p(rhs, *this);
            Fn f = p.parse_top();
            m_consts[name] = f(Env{0,0,0,0});
        }
        // Functions: "Name = expr"; each becomes a user function of (x,y,z,t).
        for (const std::string& line : funct) {
            auto eq = line.find('=');
            if (eq == std::string::npos) continue;
            std::string name = line.substr(0, eq);
            name.erase(0, name.find_first_not_of(" \t"));
            name.erase(name.find_last_not_of(" \t") + 1);
            std::string rhs = line.substr(eq + 1);
            Parser p(rhs, *this);
            Fn f = p.parse_top();          // resolves references to prior versions
            m_defs.push_back(f);
            m_funcIndex[name] = m_defs.size() - 1; // redefinition points future refs here
        }
        // Final iso expression.
        Parser p(fxyz, *this);
        m_root = p.parse_top();
    } catch (const std::exception& ex) {
        error = ex.what();
        m_root = nullptr;
        return false;
    }
    return true;
}

} // namespace Slic3r
