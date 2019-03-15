#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <functional>

struct Root
{
public:
    Root( double value, double min, double max ) : 
        value(value), min(min), max(max) 
    {
    }

    friend std::ostream& operator<<( std::ostream & os, const Root & r );
        
    double value;
    double min;
    double max;
};

std::ostream& operator<< (std::ostream & os, const Root & r)
{
    os << std::fixed << std::setprecision(13);
    os << "(value) " << r.value << "; (min) " << r.min << "; (max) " << r.max;
    return os;
} 

template<typename T>
class Functor
{
public:
    Functor( T t )
    {
        calls = 0;
        f = t;
    }

    T & operator()()
    {
        ++calls;
        return f;
    }

    int get_calls()
    {
        return calls;
    }

    void zero_out()
    {
        calls = 0;
    }

private:
    int calls;
    T f;
};


double f( const double x )
{
    return std::exp(x) - std::pow(x, 3.0);
}

Root bisection( Functor<double (*)(double)> & f, double a, double b, const double eps )
{
    double fa = f()( a );
    double fb = f()( b );

    if ( fa * fb > 0.0 )
    {
        std::cerr << "Given interval doesn't bracket zero!" << std::endl;
        exit( 1 );
    }

    double x, fx;
    while ( std::abs(a - b) > eps )
    {
        x = 0.5 * (a + b);
        fx = f()( x );

        if ( fx * fa > 0.0 )
        {
            a = x;
            fa = fx;
        }
        else
        {
            b = x;
            fb = fx;
        }
    }

    return Root( 0.5 * (a + b), a, b ); 
}

double newton_method( Functor<double(*)(double)> & f, double a, double b, const double eps )
{
    double x = 0.5 * (a + b);
    double diff, h, fx;
    
    do
    {
        h = 0.001 * ( std::abs(x) + 1.0 );
        fx = f()( x );
        diff = h * fx / ( f()(x + h) - fx );
        x -= diff;
    }
    while ( std::abs(diff) > eps );

    return x; 
}

double trisection( Functor<double(*)(double)> & f, double a, double b, const double eps )
{
    double fa = f()( a );
    double fb = f()( b );

    if ( fa * fb > 0.0 )
    {
        std::cerr << "Given interval doesn't bracket zero!" << std::endl;
        exit( 1 );
    }

    double x1, fx1;
    double x2, fx2;

    while ( b - a > eps )
    {
        if ( std::abs(fa) < std::abs(fb) )
        {
            x1 = a + (b - a) / 3.0;
            fx1 = f()( x1 );

            if ( fa * fx1 < 0.0 )
            {
                b = x1;
                fb = fx1;
            }
            else
            {
                x2 = b - (b - a) / 3.0;
                fx2 = f()( x2 );

                if ( fx1 * fx2 < 0.0 )
                {
                    a = x1;
                    fa = fx1;
                    b = x2;
                    fb = fx2;
                }
                else
                {
                    a = x2;
                    fa = fx2;
                }
            }
        }
        else
        {
            x1 = b - (b - a) / 3.0;
            fx1 = f()( x1 );
        
            if ( fb * fx1 < 0.0 )
            {
                a = x1;
                fa = fx1;
            }
            else
            {
                x2 = a + (b - a) / 3.0;
                fx2 = f()( x2 );

                if ( fx1 * fx2 < 0.0 )
                {
                    a = x2;
                    fa = fx2;
                    b = x1;
                    fb = fx1;
                }
                else
                {
                    b = x2;
                    fb = fx2;
                }
            }
        }
    }

    return 0.5 * (a + b);
}

double Interpolate2( double x1, double x2, double fx1, double fx2 )
// inverse linear interpolation 
{
    return (x1 * fx2 - x2 * fx1) / (fx2 - fx1);
}

double trisection_plus( Functor<double(*)(double)> & f, double a, double b, const double eps, const double fx_eps )
{
    double fa = f()( a );
    double fb = f()( b );

    if ( fa * fb > 0.0 )
    {
        std::cerr << "Given interval doesn't bracket zero!" << std::endl;
        exit( 1 );
    }

    double lastA, lastB;
    double x1, fx1;
    double x2, fx2;
    double x3, fx3;

    do
    {
        lastA = a;
        lastB = b;

        if ( std::abs(fa) < std::abs(fb) )
        {
            x1 = a + (b - a) / 3.0;
            fx1 = f()( x1 );

            // case 1: [a, x1] has the root
            if ( fx1 * fa < 0.0 )
            {
                x3 = Interpolate2( a, x1, fa, fx1 );
                fx3 = f()( x3 );

                if ( fa * fx3 < 0.0 )
                {
                    b = x3;
                    fb = fx3;
                }
                else
                {
                    a = x3;
                    fa = fx3;
                    b = x1;
                    fb = fx1;
                }
            }
            else
            {
                x2 = a + 2.0 * (b - a) / 3.0;
                fx2 = f()( x2 );

                // case 2: [x1, x2] has root
                if ( fx1 * fx2 < 0.0 )
                {
                    x3 = Interpolate2( x1, x2, fx1, fx2 );
                    fx3 = f()( x3 );

                    if ( fx1 * fx3 < 0.0 )
                    {
                        a = x1;
                        fa = fx1;
                        b = x3;
                        fb = fx3;
                    }
                    else
                    {
                        a = x3;
                        fa = fx3;
                        b = x2;
                        fb = fx2;
                    }
                }
                else
                {
                    // case 2: [x2, b] has the root
                    x3 = Interpolate2( x2, b, fx2, fb );
                    fx3 = f()( x3 );

                    if ( fx2 * fx3 < 0.0 )
                    {
                        a = x2;
                        fa = fx2;
                        b = x3;
                        fb = fx3;
                    }
                    else
                    {
                        a = x3;
                        fa = fx3;
                    }
                }
            }
        }
        else
        {
            x1 = a + 2.0 * (b - a) / 3.0;
            fx1 = f()( x1 );

            // case 4: [x1, b] has the root
            if ( fx1 * fb < 0.0 )
            {
                x3 = Interpolate2( x1, b, fx1, fb );
                fx3 = f()( x3 );

                if ( fx1 * fx3 < 0.0 )
                {
                    a = x1;
                    fa = fx1;
                    b = x3;
                    fb = fx3;
                }
                else
                {
                    a = x3;
                    fa = fx3;
                }
            }
            else
            {
                x2 = a + (b - a) / 3.0;
                fx2 = f()( x2 );

                // case 5: [x2, x2] has the root
                if( fx1 * fx2 < 0.0 )
                {
                    x3 = Interpolate2( x1, x2, fx1, fx2 );
                    fx3 = f()( x3 );

                    if ( fx1 * fx3 < 0.0 )
                    {
                        a = x1;
                        fa = fx1;
                        b = x3;
                        fb = fx3;
                    }
                    else
                    {
                        a = x3;
                        fa = fx3;
                        b = x2;
                        fb = fx2;
                    }
                }
                else
                {
                    // case 6: [a, x2] has the root
                    x3 = Interpolate2( a, x2, fa, fx2 );
                    fx3 = f()( x3 );

                    if ( fa * fx3 < 0.0 )
                    {
                        b = x3;
                        fb = fx3;
                    }
                    else
                    {
                        a = x3;
                        fa = x3;
                        b = x2;
                        fb = fx2;
                    }
                }
            }
        }

        if ( a > b )
        {
            std::swap(a, b);
            std::swap(fa, fb);
            std::swap(lastA, lastB);
        }

        if ( (lastA != a) && (std::abs(a - lastA) < eps) )
            break;
        if ( (lastB != b) && (std::abs(b - lastB) < eps) )
            break;
    }
    while( (std::abs(a - b) > eps) || (std::abs(fa) > fx_eps) || (std::abs(fb) > fx_eps) );

    if ( std::abs(fa) < std::abs(fb) )
        return a;
    else
        return b;
}

int main()
{
    const double a = 1.0;
    const double b = 2.0;
    const double eps = 1.0e-10;
    const double fx_eps = 1.0e-7;

    Functor<double (*)(double)> counting_f( f );

    Root r_bisection = bisection( counting_f, a, b, eps ); 
    std::cout << "Bisection root: " << r_bisection << std::endl;
    std::cout << "Number of calls: " << counting_f.get_calls() << std::endl;

    counting_f.zero_out();

    double r_newton = newton_method( counting_f, a, b, eps );
    std::cout << "Newton root: " << r_newton << std::endl;
    std::cout << "Number of calls: " << counting_f.get_calls() << std::endl;

    counting_f.zero_out();

    double r_trisection = trisection( counting_f, a, b, eps );
    std::cout << "Trisection root: " << r_trisection << std::endl;
    std::cout << "Number of calls: " << counting_f.get_calls() << std::endl;

    counting_f.zero_out();

    double r_plus = trisection_plus( counting_f, a, b, eps, fx_eps );
    std::cout << "Trisection-plus root: " << r_plus << std::endl;
    std::cout << "Number of calls: " << counting_f.get_calls() << std::endl;

    return 0;
}

