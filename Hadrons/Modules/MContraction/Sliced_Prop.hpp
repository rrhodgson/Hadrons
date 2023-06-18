/*
 * Sliced_Prop.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Lanny91 <andrew.lawson@gmail.com>
 * Author: ferben <ferben@debian.felix.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef Hadrons_MContraction_Sliced_Prop_hpp_
#define Hadrons_MContraction_Sliced_Prop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

class Sliced_PropPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Sliced_PropPar,
                                    std::string,    q,
                                    std::string,    source,
                                    std::string,    sink,
                                    int,            t,
                                    std::string,    output);
};

template <typename FImpl>
class TSliced_Prop: public Module<Sliced_PropPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string,  source,
                                        std::string,  sink,
                                        int, t);
    };
    typedef Correlator<Metadata,SitePropagator> Result;
public:
    // constructor
    TSliced_Prop(const std::string name);
    // destructor
    virtual ~TSliced_Prop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(Sliced_Prop, TSliced_Prop<FIMPL>, MContraction);

/******************************************************************************
 *                 TSliced_Prop implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TSliced_Prop<FImpl>::TSliced_Prop(const std::string name)
: Module<Sliced_PropPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TSliced_Prop<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().q};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TSliced_Prop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl>
std::vector<std::string> TSliced_Prop<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSliced_Prop<FImpl>::setup(void)
{
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TSliced_Prop<FImpl>::execute(void)
{
    // LOG(Message) << "Computing mesonic weak 3pt contractions, non-eye topologies" << std::endl;
    // LOG(Message) << "gIn : " << par().gammaIn  << std::endl;
    // LOG(Message) << "gOut: " << par().gammaOut << std::endl;
    // LOG(Message) << "q_loop  : " << par().q_loop   << std::endl;
    // LOG(Message) << "qs_left : " << par().qs_left  << std::endl;
    // LOG(Message) << "qd_right: " << par().qd_right << std::endl;
    // LOG(Message) << "qu      : " << par().qu       << std::endl;

    std::vector<Result> result;
    Result              r;
    auto                &q   = envGet(SlicedPropagator, par().q);
    // SitePropagator      q_t  = q[par().t];
    
    r.info.source  = par().source;
    r.info.sink    = par().sink;
    r.info.t       = par().t;
    
    r.corr.clear();
   
    for (int t=0; t<q.size(); t++) {
        if (par().t < 0 || t == par().t)
        r.corr.push_back( q[t] );
    }    
    result.push_back(r);

    saveResult(par().output, "Sliced_Prop", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_Sliced_Prop_hpp_
