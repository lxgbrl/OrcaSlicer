#ifndef slic3r_ImplicitExpr_hpp_
#define slic3r_ImplicitExpr_hpp_

#include <string>
#include <vector>
#include <map>
#include <functional>

namespace Slic3r {

// A small interpreter for implicit scalar fields f(x,y,z,t), compatible with a
// useful subset of MathMod's Iso3D language. Used by the custom scriptable
// infill engine (3D Volume mode, type "expr"/"mathmod") to turn user formulas
// into per-layer iso-contours.
//
// Semantics (mirrors MathMod):
//   * Every user function defined in `funct` is implicitly a function of
//     (x,y,z,t). A call Name(a,b,c,d) evaluates Name's body with x,y,z,t
//     rebound to the evaluated arguments (missing args inherit the caller's).
//   * Definitions are processed in order; a redefinition of a name binds, on
//     its right-hand side, to the *previous* version of that name.
//   * `consts` entries ("k=3/2") are evaluated once (no variables).
//   * Supported: + - * / ^, unary -, comparisons (< > <= >= == !=) yielding
//     1/0, && ||, parentheses, numbers (incl. 1/1000, 1e-3), the variables
//     x y z t, pi/PI, and builtins:
//       sin cos tan asin acos atan atan2 sinh cosh tanh sqrt cbrt abs exp
//       ln log log10 pow mod fmod sign floor ceil round min max
//       if(cond,a,b)  (lazy in the branches)
class ImplicitExpr
{
public:
    struct Env { double x, y, z, t; };
    using Fn = std::function<double(const Env&)>;

    ImplicitExpr() = default;

    // MathMod-style program: const list, funct list, and the final Fxyz expression.
    bool compile(const std::vector<std::string>& consts,
                 const std::vector<std::string>& funct,
                 const std::string&               fxyz,
                 std::string&                     error);

    // Single formula in x,y,z(,t), with optional "name=value" constants.
    bool compile_formula(const std::string&              formula,
                         const std::vector<std::string>& consts,
                         std::string&                    error);

    bool   valid() const { return bool(m_root); }
    double eval(double x, double y, double z, double t = 0.0) const
    {
        return m_root ? m_root(Env{x, y, z, t}) : 0.0;
    }

private:
    // Recursive-descent parser that compiles directly to closures.
    struct Parser;
    friend struct Parser;

    Fn                          m_root;
    std::vector<Fn>             m_defs;       // compiled user-function versions
    std::map<std::string, std::size_t> m_funcIndex; // name -> latest def index
    std::map<std::string, double>      m_consts;
};

} // namespace Slic3r

#endif // slic3r_ImplicitExpr_hpp_
