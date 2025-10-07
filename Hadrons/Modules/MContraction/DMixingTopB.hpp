/*
 * WeakEye3pt.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2023
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Felix Erben <felix.erben@cern.ch>
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
#ifndef Hadrons_MContraction_DMixingTopB_hpp_
#define Hadrons_MContraction_DMixingTopB_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopB                                        *
 *                (Fig. 4 (B) in arxiv:2504.16189)  
 *                     qCL   ┌───┐       qUR
 *                 /----<----| r |---------<---------\
 *                /          └───┘----\               \
 *               /             \       \ qInt2         \
 *           g5 *               \       \               * g5
 *               \         qInt1 \       \             /   
 *                \               \----┌───┐          /       
 *                 \--------->---------| r'|---->----/  
 *                          qUL        └───┘  qCR
 *
 * four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 * rr'=11: tr()
 * rr'=12: tr()
 * rr'=21: tr()
 * rr'=22: tr()
 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DMixingTopBPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopBPar,
                                    std::string,    qULeft,
                                    std::string,    qCLeft,
                                    std::string,    qURight,
                                    std::string,    qCRight,
                                    std::string,    qInt1,
                                    std::string,    qInt2,
                                    std::string,    output);
};

template <typename FImpl>
class TDMixingTopB: public Module<DMixingTopBPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TDMixingTopB(const std::string name);
    // destructor
    virtual ~TDMixingTopB(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
};

MODULE_REGISTER_TMP(DMixingTopB, TDMixingTopB<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopB implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopB<FImpl>::TDMixingTopB(const std::string name)
: Module<DMixingTopBPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopB<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft, 
	                           par().qCLeft,
	                           par().qURight,
	                           par().qCRight,
	                           par().qInt1,
	                           par().qInt2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopB<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopB<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopB<FImpl>::setup(void)
{
    envTmpLat(ComplexField, "corr");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopB<FImpl>::execute(void)
{

}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopB_hpp_
