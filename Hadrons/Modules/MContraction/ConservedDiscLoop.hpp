/*
 * ConservedDiscLoop.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Raoul Hodgson <raoul.hodgson@desy.de>
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
#ifndef Hadrons_MContraction_ConservedDiscLoop_hpp_
#define Hadrons_MContraction_ConservedDiscLoop_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/*
  Ward Identity contractions
 -----------------------------
 
 * options:
 - prop:       propagator. Must match the action, i.e. 5D action needs 5D propagator
 - action:     action module used for propagator solution (string)
 - source:     source module for the quark, used to remove contact terms (string)
 - mass:       mass of quark (double)
 - output:     filename for output (string)
*/

/******************************************************************************
 *                              ConservedDiscLoop                                  *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class ConservedDiscLoopPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ConservedDiscLoopPar,
                                    std::string, prop1,
                                    std::string, prop2,
                                    std::string, action,
                                    std::string, source,
                                    std::vector<int>, mus,
                                    std::string, mom,
                                    std::string, output);
};

template <typename FImpl>
class TConservedDiscLoop: public Module<ConservedDiscLoopPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Result: Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Result,
                                        unsigned int, mu,
                                        std::vector<int>, mom,
                                        std::vector<Complex>, corr);
    };

public:
    // constructor
    TConservedDiscLoop(const std::string name);
    // destructor
    virtual ~TConservedDiscLoop(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);

protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls1_;
    unsigned int Ls2_;
    std::vector<int> mom_;
};

MODULE_REGISTER_TMP(ConservedDiscLoop, TConservedDiscLoop<FIMPL>, MContraction);
MODULE_REGISTER_TMP(ZConservedDiscLoop, TConservedDiscLoop<ZFIMPL>, MContraction);

/******************************************************************************
 *                     TConservedDiscLoop implementation                           *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TConservedDiscLoop<FImpl>::TConservedDiscLoop(const std::string name)
: Module<ConservedDiscLoopPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TConservedDiscLoop<FImpl>::getInput(void)
{
    return { par().prop1, par().prop2, par().action, par().source };
}

template <typename FImpl>
std::vector<std::string> TConservedDiscLoop<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

template <typename FImpl>
std::vector<std::string> TConservedDiscLoop<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;
    
    if (!par().output.empty())
        output.push_back(resultFilename(par().output));
    
    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConservedDiscLoop<FImpl>::setup(void)
{
    // The propagator can be 4d or 5d, but must match the action
    const unsigned int ActionLs_{ env().getObjectLs(par().action) };
    Ls1_ = env().getObjectLs( par().prop1 );
    Ls2_ = env().getObjectLs( par().prop2 );
    if (Ls1_ != Ls2_ or Ls1_ != ActionLs_)
    {
        std::string sError{ "Ls mismatch: propagator Ls="};
        sError.append( std::to_string( Ls1_ ) );
        sError.append( ", action Ls=" );
        sError.append( std::to_string( ActionLs_ ) );
        HADRONS_ERROR(Size, sError);
    }

    auto parse_vector = [](const std::string &vec, int dim,
            const std::string &desc)
    {
        std::vector<int> res = strToVec<int>(vec);
        if(res.size() != dim) {
            HADRONS_ERROR(Size, desc + " has "
                    + std::to_string(res.size()) + " instead of "
                    + std::to_string(dim) + " components");
        }
        return res;
    };
    mom_      = parse_vector(par().mom, env().getNd()-1, "momentum");

    envTmpLat(PropagatorField, "tmp");
    envTmpLat(ComplexField, "tmp_scalar");
    envTmpLat(LatticeComplex, "coor");
    envTmpLat(LatticeComplex, "ph");
    envTmpLat(LatticeComplex, "c");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TConservedDiscLoop<FImpl>::execute(void)
{
    // LOG(Message) << "Performing Ward Identity checks for propagator " << par().prop << std::endl;
    auto &prop1 = envGet(PropagatorField, par().prop1);
    auto &prop2 = envGet(PropagatorField, par().prop2);
    LOG(Message) << "Action " << par().action << std::endl;
    auto &act = envGet(FMat, par().action);
    LOG(Message) << "Physical source " << par().source << std::endl;
    auto &phys_source = envGet(PropagatorField, par().source);
    // Gamma g5(Gamma::Algebra::Gamma5);
    // Gamma gT(Gamma::Algebra::GammaT);


    envGetTmp(LatticeComplex, ph);
    ph = Zero();
    if (mom_[0] != 0 || mom_[1] != 0 || mom_[2] != 0) 
    {
        LOG(Message) << "Adding momentum phase " << mom_ << std::endl;

        envGetTmp(LatticeComplex, coor);
        for(unsigned int mu = 0; mu < 3; mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (1.*mom_[mu]/env().getDim(mu))*coor;
        }
    }
    Complex i(0.0,1.0);
    ph = exp((Real)(2*M_PI)*i*ph);
    

    // Create results = zero
    std::vector<Result> result;
    result.resize(par().mus.size());

    envGetTmp(PropagatorField, tmp);
    envGetTmp(ComplexField, tmp_scalar);
    std::vector<TComplex> buf;

    for (int i=0; i<par().mus.size(); i++)
    {
        int mu = par().mus[i];

        act.ContractConservedCurrent(prop1, prop2, tmp, phys_source, Current::Vector, mu);
        tmp_scalar = trace(tmp)*ph;

        sliceSum(tmp_scalar, buf, Tp);
        result[i].corr.resize(buf.size());
        for (unsigned int t = 0; t < buf.size(); ++t)
            result[i].corr[t] = TensorRemove(buf[t]);

        result[i].mu = mu;
        result[i].mom = mom_;
    }

    saveResult(par().output, "ConservedDiscLoop", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_ConservedDiscLoop_hpp_
