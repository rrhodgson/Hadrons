/*
 * ScaledDWF.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef Hadrons_MAction_ScaledDWF_hpp_
#define Hadrons_MAction_ScaledDWF_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Modules/MAction/FermionActionModule.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                      Scaled domain wall fermion                            *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MAction)

class ScaledDWFPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(ScaledDWFPar,
                                    std::string , gauge,
                                    unsigned int, Ls,
                                    double      , mass,
                                    double      , M5,
                                    double      , scale,
                                    std::string , boundary,
                                    std::string , twist);
};

template <typename FImpl>
class TScaledDWF: public FermionActionModule<ScaledDWFPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TScaledDWF(const std::string name);
    // destructor
    virtual ~TScaledDWF(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(ScaledDWF, TScaledDWF<FIMPL>, MAction);
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
MODULE_REGISTER_TMP(ScaledDWFF, TScaledDWF<FIMPLF>, MAction);
#endif

/******************************************************************************
 *                      TScaledDWF implementation                             *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TScaledDWF<FImpl>::TScaledDWF(const std::string name)
: FermionActionModule<ScaledDWFPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TScaledDWF<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().gauge};

    if ((!isVector<Real>(par().twist)) && (!par().twist.empty()))
    {
        in.push_back(par().twist);
    }
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TScaledDWF<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TScaledDWF<FImpl>::setup(void)
{
    LOG(Message) << "Setting up scaled domain wall fermion matrix with m= "
                 << par().mass << ", M5= " << par().M5 << ", Ls= " << par().Ls 
                 << ", scale= " << par().scale
                 << " using gauge field '" << par().gauge << "'"
                 << std::endl;

    auto &U    = envGet(GaugeField, par().gauge);
    auto &g4   = *envGetGrid(FermionField);
    auto &grb4 = *envGetRbGrid(FermionField);
    auto &g5   = *envGetGrid(FermionField, par().Ls);
    auto &grb5 = *envGetRbGrid(FermionField, par().Ls);
    typename ScaledShamirFermion<FImpl>::ImplParams implParams;
    parseBoundary(implParams);
    envCreateDerived(FMat, ScaledShamirFermion<FImpl>, getName(), par().Ls, U, g5,
                     grb5, g4, grb4, par().mass, par().M5, par().scale,
                     implParams);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TScaledDWF<FImpl>::execute(void)
{}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MAction_ScaledDWF_hpp_
