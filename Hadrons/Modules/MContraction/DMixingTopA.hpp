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
#ifndef Hadrons_MContraction_DMixingTopA_hpp_
#define Hadrons_MContraction_DMixingTopA_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopA                                        *
 *                (Fig. 4 (A) in arxiv:2504.16189)  
 *                   qCL        qInt1        qUR
 *                 /--<--\   /----<----\   /--<--\
 *                /       \ /           \ /       \
 *               /       ┌───┐         ┌───┐       \
 *           g5 *        | r |         | r'|        * g5
 *               \       └───┘         └───┘       /
 *                \       / \           / \       /
 *                 \-->--/   \---->----/   \-->--/
 *                   qUL        qInt1        qCR
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

class DMixingTopAPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopAPar,
                                    std::string,    qULeft,
                                    std::string,    qCLeft,
                                    std::string,    qURight,
                                    std::string,    qCRight,
                                    std::string,    qInt1,
                                    std::string,    qInt2,
                                    std::string,    output);
};

template <typename FImpl>
class TDMixingTopA: public Module<DMixingTopAPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string, rr,
                                        std::string, parity,
                                        std::string, eta_max);
    };
    typedef Correlator<Metadata, std::vector<Complex>> Result;
public:
    // constructor
    TDMixingTopA(const std::string name);
    // destructor
    virtual ~TDMixingTopA(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual PropagatorField GH_VVAA_cap(const PropagatorField &prop, int r);
};

MODULE_REGISTER_TMP(DMixingTopA, TDMixingTopA<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopA implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopA<FImpl>::TDMixingTopA(const std::string name)
: Module<DMixingTopAPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopA<FImpl>::getInput(void)
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
std::vector<std::string> TDMixingTopA<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopA<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

template <typename FImpl>
typename TDMixingTopA<FImpl>::PropagatorField TDMixingTopA<FImpl>::GH_VVAA_cap(const TDMixingTopA<FImpl>::PropagatorField &prop, int r)
{
    assert(r==1 or r==2);

    GridBase *grid = envGetGrid(FermionField);

    std::array<Gamma, 8> GHs{Gamma(Gamma::Algebra::GammaX),
                             Gamma(Gamma::Algebra::GammaY),
                             Gamma(Gamma::Algebra::GammaZ),
                             Gamma(Gamma::Algebra::GammaT),
                             Gamma(Gamma::Algebra::GammaXGamma5),
                             Gamma(Gamma::Algebra::GammaYGamma5),
                             Gamma(Gamma::Algebra::GammaZGamma5),
                             Gamma(Gamma::Algebra::GammaTGamma5)};

    SitePropagator spId(1.0);

    PropagatorField GPropG_VVAA(grid);
    GPropG_VVAA = Zero();
    for (int g = 0; g < GHs.size(); g++)
    {
        Gamma GH = GHs[g];
        if (r == 1)
        {
            GPropG_VVAA += spId * GH * trace(prop * GH);
        }
        else
        {
            GPropG_VVAA += GH * prop * GH;
        }
    }
    return GPropG_VVAA;
};

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopA<FImpl>::setup(void)
{
    envTmpLat(PropagatorField, "qcul");
    envTmpLat(PropagatorField, "qcur");

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopA<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology D" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    LOG(Message) << "qInt1   : " << par().qInt1 << std::endl;
    LOG(Message) << "qInt2   : " << par().qInt2 << std::endl;
      
    std::vector<Result> result;
    Result res;
    
    const int Nt{env().getDim(Tdir)};
    GridCartesian *grid = envGetGrid(FermionField);
    
    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);
    auto &qi1 = envGet(std::vector<PropagatorField *>, par().qInt1);
    auto &qi2 = envGet(std::vector<PropagatorField *>, par().qInt2);

    envGetTmp(PropagatorField, qcul);
    envGetTmp(PropagatorField, qcur);

    Gamma g5(Gamma::Algebra::Gamma5);   

    qcul = qcl * g5 * g5 * adj(qul) * g5;
    qcur = qcr * g5 * g5 * adj(qur) * g5;

}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopA_hpp_
