/*
 * PropGamma.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_PropGamma_hpp_
#define Hadrons_MContraction_PropGamma_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

BEGIN_MODULE_NAMESPACE(MContraction)

class PropGammaPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(PropGammaPar,
                                    std::string,    src,
                                    std::string,    snk,
                                    std::vector<std::string>, moms,
                                    std::string, gamma,
                                    std::string,    output);
};

template <typename FImpl>
class TPropGamma: public Module<PropGammaPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string,    mom,
                                        Gamma::Algebra, gamma);
    };
    typedef Correlator<Metadata,SitePropagator> Result;
public:
    // constructor
    TPropGamma(const std::string name);
    // destructor
    virtual ~TPropGamma(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);

    void parseGammaString(std::vector<Gamma::Algebra> &gammaList);
};

MODULE_REGISTER_TMP(PropGamma, TPropGamma<FIMPL>, MContraction);

/******************************************************************************
 *                 TPropGamma implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TPropGamma<FImpl>::TPropGamma(const std::string name)
: Module<PropGammaPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TPropGamma<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().src,par().snk};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TPropGamma<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {};
    
    return out;
}

template <typename FImpl>
std::vector<std::string> TPropGamma<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output = {resultFilename(par().output)};
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropGamma<FImpl>::setup(void)
{
}

template <typename FImpl>
void TPropGamma<FImpl>::parseGammaString(std::vector<Gamma::Algebra> &gammaList)
{
    gammaList.clear();
    // Determine gamma matrices to insert at source/sink.
    if (par().gamma.compare("all") == 0)
    {
        // Do all contractions.
        for (unsigned int i = 1; i < Gamma::nGamma; i += 2)
        {
            gammaList.push_back((Gamma::Algebra)i);
        }
    }
    else
    {
        // Parse individual contractions from input string.
        gammaList = strToVec<Gamma::Algebra>(par().gamma);
    } 
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TPropGamma<FImpl>::execute(void)
{
    // LOG(Message) << "Computing mesonic weak 3pt contractions, non-eye topologies" << std::endl;
    // LOG(Message) << "gIn : " << par().gammaIn  << std::endl;
    // LOG(Message) << "gOut: " << par().gammaOut << std::endl;
    // LOG(Message) << "q_loop  : " << par().q_loop   << std::endl;
    // LOG(Message) << "qs_left : " << par().qs_left  << std::endl;
    // LOG(Message) << "qd_right: " << par().qd_right << std::endl;
    // LOG(Message) << "qu      : " << par().qu       << std::endl;


    unsigned int nd    = env().getDim().size();

    std::vector<Result> result;
    Result              r;

    auto                &src   = envGet(PropagatorField, par().src);
    auto                &snk   = envGet(PropagatorField, par().snk);

    PropagatorField     tmp(envGetGrid(PropagatorField));
    SlicedPropagator    res;

    LatticeComplex coor(envGetGrid(LatticeComplex));
    LatticeComplex ph(envGetGrid(LatticeComplex));

    Gamma g5(Gamma::Algebra::Gamma5);
    std::vector<Gamma::Algebra> gammaList;
    parseGammaString(gammaList);

    Complex i(0.0,1.0);

    for (const auto& mom : par().moms) {
        std::vector<Real> p = strToVec<Real>(mom);
        if (p.size() != nd-1) {
            HADRONS_ERROR(Size, "momentum number of components different from " + std::to_string(nd-1));
        }

        ph = Zero();
        for(unsigned int mu = 0; mu < p.size(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);

        for (const auto& gamma : gammaList) {
            Gamma G(gamma);

            tmp = g5*adj(snk)*g5*G*src;
            tmp *= ph;

            sliceSum(tmp, res, Tp);
            r.corr = res;
            r.info.gamma = gamma;
            r.info.mom = mom;
            result.push_back(r);
        }
    }

    saveResult(par().output, "PropGamma", result);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_PropGamma_hpp_
