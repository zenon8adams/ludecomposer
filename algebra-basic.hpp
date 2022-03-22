#ifndef __ALGEBRA_BASIC_HPP
#define __ALGEBRA_BASIC_HPP

#include <vector>
#include <set>
#include <string>
#include <initializer_list>
#include <iterator>
#include <unordered_map>
#include <numeric>
#include <cmath>
#include <iostream>

#include "macros.hpp"

class Variable
{

public:
    std::string name;
    int index;

    Variable( std::string nm, int idx = 1)
        : name( nm), index( idx)
    {
    }

    std::string hashVal() const
    {
        std::string result( index != 0 ? name : "");
        if( index == 1 || index == 0 || name.empty())
            return result;

        result.insert( result.end(), '^');
        result.insert( result.size(), std::to_string( index));
        
        return result;
    }

    bool operator==( const Variable& rhs) const
    {
        return hashVal() == rhs.hashVal();
    }

    bool operator<( const Variable& rhs) const
    {
        return hashVal() < rhs.hashVal(); 
    }

    value_t eval( value_t val) const
    {
#ifdef __MP
        value_t result( val.get_mpf_t());
        mpf_pow_ui( result.get_mpf_t(), result.get_mpf_t(), index);
        return result;
#else
        return std::pow( val, index);
#endif
    }

};
class Monomial
{
public:

template <typename Container = std::set<Variable>>
    explicit Monomial( value_t co = 1, const Container &tms = {})
        : coeff( co)
    {
        if( !ZERO( coeff))
            std::copy( tms.cbegin(), tms.cend(), std::inserter( term, term.end()));
    }

    Monomial( const Monomial& m)
        : coeff( m.coeff)
    {
        assign( m, *this);
    }

    Monomial& operator=( const Monomial& m)
    {
        if( this == &m)
            return *this;

        assign( m, *this);
        return *this;
    }
    
    Monomial& operator-()
    {
        coeff = -coeff;
        return *this;
    }

    auto differential( const std::string& val) const
    {
        auto diffFn = []( auto& co, auto& var, auto& result)
        {
            if( var.name.empty())
                return Monomial{ 0};

            co *= var.index--;
            if( var.index != 0)
                result.push_back( var);

            return Monomial( co, result);
        };

        return apply( val, diffFn);
    }
        
    auto integral( const std::string& val) const
    {
        auto integFn = [ &val]( auto& co, auto& var, auto& result)
        {
            if( var.name.empty())
                result.push_back( Variable( val));
            else
            {
                co /= ++var.index;
                if( var.index != 0)
                    result.push_back( var);
            }

            return Monomial( co, result);
        };

        return apply( val, integFn);
    }

    auto dump() const
    {
        return std::vector<Variable>( term.cbegin(), term.cend());
    }

    Monomial operator*( const Monomial& rhs) const
    {
        Monomial result( coeff, {});
        result.coeff *= rhs.coeff;
        if( ZERO( result.coeff))
            return result;
        std::unordered_map<std::string, double> power;
        for( const auto& unit : term)
            power[ unit.name] += unit.index;
        for( const auto& unit : rhs.term)
            power[ unit.name] += unit.index;
        std::for_each( power.cbegin(), power.cend(), [&]( auto& elem)
                { result.term.emplace( elem.first, elem.second); });
        return result;
    }

    Monomial& operator*=( const Monomial& rhs)
    {
        return *this = *this * rhs;
    }

    std::string hashVal() const
    {
        return std::accumulate( term.cbegin(), term.cend(), std::string(),
                []( const auto& l, const auto& r){ return l + r.hashVal(); });
    }

    std::string toString() const
    {
        if( EQUAL( coeff, 1))
        {
            std::string ret( hashVal());
            return ret.empty() ? "1" : ret;
        }
#ifdef __MP
        mpf_t copy;
        mpf_init_set( copy, coeff.get_mpf_t());
        char *str = NULL;
        gmp_asprintf( &str, "%Fg", copy);
        std::string out( str);
        free( str);
        mpf_clear( copy);
#endif
        return
#ifdef __MP
         
        out
#else
        std::to_string( coeff)
#endif
        + hashVal();
    }

    void setCoeff( value_t new_co)
    {
        coeff = new_co;
    }

    value_t coefficient() const
    {
        return coeff;
    }

    size_t unknown() const
    {
        return term.size();
    }

    Monomial eval( const result_t& range) const
    {
        auto cf = coeff;
        std::set<Variable> new_t;
        for( const auto& elem : term)
        {
            auto val = range.find( elem.name);
            if( val == range.cend())
                continue;
            if( !val->second.second)
            {
                new_t.insert( elem);
                continue;
            }
            
            cf *= elem.eval( val->second.first);
        }

        return Monomial( cf, new_t);
    }

private:

    template <typename Function>
    Monomial apply( const std::string& val, Function fun) const
    {
        std::vector<Variable> result;
        value_t co = coeff;
        Variable indep( {});
        std::copy_if( term.cbegin(), term.cend(),
                std::back_inserter( result),
                [ &]( const auto& var)
                {
                   if( var.name == val)
                   {
                        indep = var;
                        return false;
                    }
                    
                    return true;
                } );

        return fun( co, indep, result);
    }

    static void assign( const Monomial& src, Monomial& dest)
    {
        if( ZERO( src.coeff))
        {
            dest.term.clear();
            return;
        }

        dest.coeff = src.coeff;
        std::copy( src.term.cbegin(), src.term.cend(),
                std::inserter( dest.term, dest.term.end()));
    }

    value_t coeff;
    std::set<Variable> term;
};

class AlUnit
{
public:

    AlUnit( const Monomial& mon)
    {
        if( !ZERO( mon.coefficient()))
            combo.push_back( mon);
    }

    AlUnit( const std::initializer_list<Monomial>& mns = {})
    {
        assign( mns, *this);
    }

    AlUnit( const AlUnit& dupl)
    {
        assign( dupl.combo, *this);
    }

    AlUnit& operator=( const AlUnit& dupl)
    {
        if( this == &dupl)
            return *this;

        assign( dupl.combo, *this);
        return *this;
    }

    AlUnit operator+( const AlUnit& rhs) const
    {
        AlUnit result;
        std::string hv;
        std::unordered_map<std::string, value_t> sum;
        std::unordered_map<std::string, Monomial> cache;

        for( const auto& term : combo)
        {
            sum[ hv = term.hashVal()] = term.coefficient();
            cache[ hv] = term;
        }

        for( const auto& term : rhs.combo)
        {
            sum[ hv = term.hashVal()] += term.coefficient();
            cache[ hv] = term;
        }

        decltype(sum)::const_iterator sum_s = sum.cbegin(),
                                      sum_e = sum.cend();

        while( sum_s != sum_e)
        {
            cache[ sum_s->first].setCoeff( sum_s->second);
            result.combo.push_back( cache[ sum_s->first]);
            ++sum_s;
        }
        return AlUnit(result);
    } 

    AlUnit& operator-()
    {
        for( auto& term : combo)
            term = -term;
        return *this;
    }

    AlUnit operator-( AlUnit rhs)
    {
        return *this + -rhs;
    }

    AlUnit& operator-=( const AlUnit& rhs)
    {
        return *this = *this - rhs;
    }

    AlUnit& operator+=( const AlUnit& rhs)
    {
        return *this = *this + rhs;
    }

    friend AlUnit operator*( value_t scale, const AlUnit& rhs)
    {
        return AlUnit(Monomial( scale)) * rhs;
    }

    AlUnit operator*( const AlUnit& rhs) const
    {
        AlUnit result;
        for( const auto& lt : combo)
            for( const auto& rt : rhs.combo)
                result += AlUnit( lt * rt);

        return AlUnit(result);
    }

    AlUnit operator*( value_t scale)
    {
        return scale * *this;
    }

    AlUnit& operator*=( const AlUnit& rhs)
    {
        return *this = *this * rhs;
    }

    std::string hashVal() const
    {
        return std::accumulate( combo.cbegin(), combo.cend(), std::string(),
                []( const auto& l, const auto& r){ return l + r.hashVal(); });
    }

    //TODO: Complete single variable Maclaurin series.

    AlUnit sVarMaclaurinSeries( const std::string& indep, value_t idx, size_t iteration) const
    {
        AlUnit current{ *this}, result;
        result_t values{ { indep, { 0, true}}};
        value_t fact = 1, prefix = 1;

        for( size_t i = 0; i < iteration; ++i)
        {
            std::cout << "Cur: " << current.toString()
                      << ", Prefix: " << prefix
                      << ", Pow: " << idx <<'\n';
            current *= current.differential( indep);
            prefix  *= idx;
            idx -= 1;
            result += prefix * current.eval( values).value()
                                    *
                            fact    *    AlUnit{ Monomial( 1, { Variable( indep, i)})};
            fact *=(int)( 1 / ( i + 1));
        }

        return result;
    }

    AlUnit differential( const std::string& indep) const
    {
        AlUnit result;

        std::transform( combo.cbegin(), combo.cend(),
                std::back_inserter( result.combo),
                [ &indep]( const auto& term)
                {
                    return term.differential( indep);
                } );

        return AlUnit( result);
    }

    AlUnit integral( const std::string& indep) const
    {
        AlUnit result;

        std::transform( combo.cbegin(), combo.cend(),
                std::back_inserter( result.combo),
                [ &indep]( const auto& term)
                {
                    return term.integral( indep);
                } );

        return AlUnit( result);
    }

    AlUnit eval( const result_t& soln)
    {
        AlUnit result;
        for( const auto& term : combo)
            result += term.eval( soln);

        return AlUnit(result);
    }

    value_t value() const
    {
        if( !unknown().empty())
            throw std::logic_error( "Value not defined!");
        if( combo.empty())
            return 0.0;

        return combo.front().coefficient();
    }

    std::vector<Monomial> unknown() const
    {
        std::vector<Monomial> result;
        for( const auto& term : combo)
            if( term.unknown() == 1)
                result.push_back( term);
        return result;
    }

    auto dump() const
    {
        return combo;
    }

    std::string toString() const
    {
        std::string result;
        decltype(combo)::const_iterator combo_s( combo.cbegin()),
                                        combo_e( combo.cend());
        if( combo_s != combo_e)
            result.insert( result.size(), (combo_s++)->toString());
        
        while( combo_s != combo_e)
        {
            if( combo_s->coefficient() > 0)
                result.insert( result.end(), '+');
            result.insert( result.size(), combo_s->toString());
            ++combo_s;
        }

        return result.empty() ? "0"s : result;
    }

private:
    template <typename Source>
    static void assign( const Source& src, AlUnit& dest)
    {
        dest.combo.clear();
        std::copy_if( src.begin(), src.end(), 
                std::back_inserter( dest.combo),
                []( const auto& elem){ return !ZERO( elem.coefficient()); });
    }

    std::vector<Monomial> combo;
};

#endif
