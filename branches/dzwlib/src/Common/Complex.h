/**   
* @file Complex.h
* @brief define the complex form of vector field
* @author dangzw
* @date Oct 20, 2011 4:06:58 PM 
* @version V1.0   
*/

#ifndef COMPLEX_H_
#define COMPLEX_H_

#include <cmath>

namespace dzw{namespace common{
/**
 * This class express 2D direction field with (cosx, sinx).
 * We can get the angle form of the vector field.
 * And this class provide some operators.
 */
class Complex
{
public:
  /**Constructor.
  * @param[in] a cosx of the direction field.
  * @param[in] b sinx of the direction field.
  * @return
  */
  Complex(double a=0, double b=0) : cosx(a), sinx(b) { }
  /**Get cosx of the direction field.
  * @return cosx of the direction field.
  */
  double cosX() const { return cosx ; }
  /**Get sinx of the direction field.
  * @return sinx of the direction field.
  */
  double sinX() const { return sinx ; }

  /** Get squared norm of the direction field.
  * @return squared norm of the direction field.
  */
  double squared_modulus() const {
    return cosx*cosx + sinx*sinx ;
  }
  /** Get norm of the direction field.
  * @return norm of the direction field.
  */
  double modulus() const {
    return ::sqrt(squared_modulus()) ;
  }
  /**Get angle of the direction field.
  * @return angle form
  */
  double angle() const {
    double result = 0.0 ;
    double m = modulus() ;
    if(m > 1e-30) {
      result = acos(cosx / modulus() );
      if (sinx < 0) {
        result = - result ;
      }
    }
    return result;
  }

  /** Normalize the direction field.
  */
  void normalize() {
    double n = modulus() ;
    double s = n < 1e-20 ? 0.0 : 1.0 / n ;
    cosx *= s ;
    sinx *= s ;
  }

  /**Get normalized direction field.
  * @return normalized direction field.
  */
  Complex normalized() const {
    Complex result = *this ;
    result.normalize() ;
    return result ;
  }
  /** operator +=
  * A += B
  * @param rhs the complex B.
  * @return
  */
  Complex& operator+=(const Complex& rhs) {
    cosx += rhs.cosx ;
    sinx += rhs.sinx ;
    return *this ;
  }
  /** operator -=
  * A -= B
  * @param rhs the complex B.
  * @return
  */
  Complex& operator-=(const Complex& rhs) {
    cosx -= rhs.cosx ;
    sinx -= rhs.sinx ;
    return *this ;
  }

private:
  double cosx ; /**< cosx of the direction field, the cos of angle x, which is between the direction field and reference frame. */
  double sinx ; /**< sinx of the direction field, as before. */
} ;

/**Operator *
 * A = a * B, A is complex, a is double type, B is also complex.
 * @param a double number.
 * @param z Complex
 * @return Result is Complex.
 */
inline Complex operator*(double a, const Complex& z) {
  return Complex(a*z.cosX(), a*z.sinX()) ;
}
/**Operator +
 * A = B + C, A, B, C are all complex type.
 * @param z1 Complex.
 * @param z2 Complex
 * @return Result is Complex.
 */
inline Complex operator+(const Complex& z1, const Complex& z2) {
  return Complex(
      z1.cosX() + z2.cosX(),
      z1.sinX() + z2.sinX()
  ) ;
}
/**Operator -
 * A = B - C, A, B, C are all complex type.
 * @param z1 Complex.
 * @param z2 Complex
 * @return Result is Complex.
 */
inline Complex operator-(const Complex& z1, const Complex& z2) {
  return Complex(
      z1.cosX() - z2.cosX(),
      z1.sinX() - z2.sinX()
  ) ;
}
/**Operator *
 * A = B * C, A, B, C are all complex type.
 * @param z1 Complex.
 * @param z2 Complex
 * @return Result is Complex.
 */
inline Complex operator*(const Complex& z1, const Complex& z2) {
  return Complex(
      z1.cosX()*z2.cosX() - z1.sinX()*z2.sinX(),
      z1.cosX()*z2.sinX() + z1.sinX()*z2.cosX()
  ) ;
}
/**Operator /
 * A = B / C, A, B, C are all complex type.
 * @param z1 Complex.
 * @param z2 Complex
 * @return Result is Complex.
 */
inline Complex operator/(const Complex& z1, const Complex& z2) {
  double d = z2.squared_modulus() ;
  return Complex(
      (z1.cosX()*z2.cosX()+z1.sinX()*z2.sinX()) / d,
      (z1.sinX()*z2.cosX() - z1.cosX()*z2.sinX()) / d
  ) ;
}
/**Operator sqrt
 * A = sqrt(B), A, B are all complex type.
 * @param z Complex.
 * @return Result is Complex.
 */
inline Complex complex_sqrt(const Complex& z){
  double r = ::sqrt(z.modulus());
  double angle = 0.5*z.angle() ;
  return Complex(r*::cos(angle),r*::sin(angle));
}
/**output operator <<
 * @param out ostream
 * @param z Complex
 * @return
 */
inline std::ostream& operator<<(std::ostream& out, const Complex& z) {
  return out << z.cosX() << " " << z.sinX() ;
}
/**output operator >>
 * @param in istream
 * @param z Complex
 * @return
 */
inline std::istream& operator>>(std::istream& in, Complex& z) {
  double cosX,imag ;
  in >> cosX >> imag ;
  z = Complex(cosX,imag) ;
  return in ;
}
}}
#endif /* COMPLEX_H_ */
