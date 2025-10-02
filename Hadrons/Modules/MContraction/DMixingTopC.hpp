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
#ifndef Hadrons_MContraction_DMixingTopC_hpp_
#define Hadrons_MContraction_DMixingTopC_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopC                                        *
 *                  (Fig. 4 (C) in arxiv:2504.16189)  
 *                 qCL                             qUR
 *               /--<--\    qLoop1                /--<--\
 *              /       \    /--\      /--\      /       \
 *             /       ┌───┐/    \    /    \┌───┐         \
 *         g5 *        | r |      |   |     | r'|          * g5
 *             \       └───┘\    /    \    /└───┘         /
 *              \       /    \--/      \--/      \       /
 *               \-->--/              qLoop2      \-->--/
 *                 qUL                              qCR
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

class DMixingTopCPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopCPar,
                                    std::string,    qULeft,
                                    std::string,    qCLeft,
                                    //std::string,    qURight,
                                    //std::string,    qCRight,
                                    std::string,    qLoop1,
                                    //std::string,    qLoop2,
                                    std::string,    output);
};

template <typename FImpl>
class TDMixingTopC: public Module<DMixingTopCPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string, r,
                                        std::string, parity,
                                        std::string, eta);
    };
    typedef Correlator<Metadata, Complex> Result;
public:
    // constructor
    TDMixingTopC(const std::string name);
    // destructor
    virtual ~TDMixingTopC(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual std::vector<Complex> contract_C_half(const LatticePropagator &prop_c, const LatticePropagator &prop_u, const LatticePropagator &loop);
    virtual std::pair<LatticePropagator, LatticePropagator> GH_VVAA_cap(const LatticePropagator &prop);
};

MODULE_REGISTER_TMP(DMixingTopC, TDMixingTopC<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopC implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopC<FImpl>::TDMixingTopC(const std::string name)
: Module<DMixingTopCPar>(name)
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopC<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft, 
                               par().qCLeft,
                               //par().qURight,
                               //par().qCRight,
                               par().qLoop1};//,
                               //par().qLoop2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopC<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopC<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

template <typename FImpl>
std::vector<Complex> TDMixingTopC<FImpl>::contract_C_half(const LatticePropagator &prop_c, const LatticePropagator &prop_u, const LatticePropagator &loop)
{
    Gamma G5(Gamma::Algebra::Gamma5);
    Gamma GT(Gamma::Algebra::GammaT);

    Gamma Gsrc = G5;

    LatticeComplex tmp = trace( prop_c * Gsrc * G5*adj(prop_u)*G5 * loop );
    std::vector<LatticeComplex::scalar_object> ret;
    sliceSum(tmp, ret, Tp);
    std::vector<Complex> ret2(ret.size());
    for (unsigned int t = 0; t < env().getDim(Tdir); ++t)
    {
        ret2[t] = TensorRemove(ret[t]);
    }
    return ret2;
};

template <typename FImpl>
std::pair<LatticePropagator, LatticePropagator> TDMixingTopC<FImpl>::GH_VVAA_cap(const LatticePropagator &prop)
{
    GridBase *grid = prop.Grid();

    std::array<Gamma, 8> GHs{Gamma(Gamma::Algebra::GammaX),
                             Gamma(Gamma::Algebra::GammaY),
                             Gamma(Gamma::Algebra::GammaZ),
                             Gamma(Gamma::Algebra::GammaT),
                             Gamma(Gamma::Algebra::GammaXGamma5),
                             Gamma(Gamma::Algebra::GammaYGamma5),
                             Gamma(Gamma::Algebra::GammaZGamma5),
                             Gamma(Gamma::Algebra::GammaTGamma5)};

    SpinColourMatrix spId = Zero();
    for (int s = 0; s < 4; s++)
    {
        for (int c = 0; c < 3; c++)
        {
            spId()(s, s)(c, c) = 1.;
        }
    }

    LatticePropagator GTrPropG_VVAA(grid);
    GTrPropG_VVAA = Zero();
    LatticePropagator GPropG_VVAA(grid);
    GPropG_VVAA = Zero();
    for (int g = 0; g < GHs.size(); g++)
    {
        Gamma GH = GHs[g];
        GTrPropG_VVAA += spId * GH * trace(prop * GH);
        GPropG_VVAA += GH * prop * GH;
    }
    return std::make_pair(GTrPropG_VVAA, GPropG_VVAA);
};

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopC<FImpl>::setup(void)
{
    GridCartesian *grid = envGetGrid(FermionField);
    envTmp(std::vector<LatticePropagator>, "GdsG_pp", 1, 2, LatticePropagator(env().getGrid()));   
    envTmpLat(ComplexField, "corr");
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopC<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology C" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qLoop1  : " << par().qLoop1 << std::endl;

    std::vector<Result> result;
    Result res;

    const int Nt{env().getDim(Tdir)};
    GridCartesian *grid = envGetGrid(FermionField);

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &ql1 = envGet(std::vector<PropagatorField *>, par().qLoop1);

    int Neta = ql1.size();

    // parity +, parity -
    std::vector<Gamma> parityG = {Gamma(Gamma::Algebra::Identity), Gamma(Gamma::Algebra::Gamma5)};
    envGetTmp(std::vector<LatticePropagator>, GdsG_pp);
    std::vector<Complex> buf(Nt);
    for (int p = 0; p < 2; p++)
    {
        res.info.parity = (p == 0) ? "+" : "-";
        for (int i = 0; i < Neta; i++)
        {
            // here one has to add ql2 if one wants them to be allowed to be different
            auto tmp = GH_VVAA_cap(*ql1[i]);
            GdsG_pp[0] = tmp.first;  // r1
            GdsG_pp[1] = tmp.second; // r2
            for (int r = 0; r < 2; r++)
            {
                res.info.r = std::to_string(r+1);
                res.info.eta = i;
                buf = contract_C_half(qcl, qul, GdsG_pp[r] * parityG[p]);
                res.corr.clear();
                res.corr = buf;
                result.push_back(res);
            }
        }
    }


}


END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopC_hpp_
