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
#ifndef Hadrons_MContraction_DMixingTopD_hpp_
#define Hadrons_MContraction_DMixingTopD_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopD                                        *
 *                    (Fig. 4 (D) in arxiv:2504.16189)
 *             qCL   ┌───┐                 qUR
 *         /----<----| r |-----------------<--------------\
 *        /          └───┘                                 \
 *       /          /     \          qLoop2 /->-\           \
 *   g5 *           \     /                /     \           * g5
 *       \           \->-/ qLoop1          \     /          /
 *        \                                 ┌───┐          /
 *         \-------------->-----------------| r'|---->----/
 *                       qUL                └───┘  qCR
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

class DMixingTopDPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopDPar,
                                    std::string, qULeft,
                                    std::string, qCLeft,
                                    std::string, qURight,
                                    std::string, qCRight,
                                    std::string, qLoop1,
                                    std::string, qLoop2,
                                    std::string, output);
};

template <typename FImpl>
class TDMixingTopD : public Module<DMixingTopDPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
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
    TDMixingTopD(const std::string name);
    // destructor
    virtual ~TDMixingTopD(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
    virtual std::vector<std::string> getOutputFiles(void);
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
    // bespoke subcontractions
    virtual SlicedPropagator contract_D_half(const PropagatorField &prop_c, const PropagatorField &prop_u, const PropagatorField &loop);
    virtual std::vector<std::vector<Complex>> contract_D(const SlicedPropagator &half_if, const SlicedPropagator &half_fi);
    virtual std::pair<PropagatorField, PropagatorField> GH_VVAA_cap(const PropagatorField &prop);
};

MODULE_REGISTER_TMP(DMixingTopD, TDMixingTopD<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopD implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopD<FImpl>::TDMixingTopD(const std::string name)
    : Module<DMixingTopDPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft,
                                   par().qCLeft,
                                   par().qURight,
                                   par().qCRight,
                                   par().qLoop1,
                                   par().qLoop2};

    return in;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};

    return out;
}

template <typename FImpl>
std::vector<std::string> TDMixingTopD<FImpl>::getOutputFiles(void)
{
    std::vector<std::string> output;

    if (!par().output.empty())
        output.push_back(resultFilename(par().output));

    return output;
}

template <typename FImpl>
typename TDMixingTopD<FImpl>::SlicedPropagator TDMixingTopD<FImpl>::contract_D_half(const TDMixingTopD<FImpl>::PropagatorField &prop_c, const TDMixingTopD<FImpl>::PropagatorField &prop_u, const TDMixingTopD<FImpl>::PropagatorField &loop)
{
    Gamma G5(Gamma::Algebra::Gamma5);

    PropagatorField tmp = G5 * adj(prop_u) * G5 * loop * prop_c;
    SlicedPropagator ret;
    sliceSum(tmp, ret, Tp);
    return ret;
};

template <typename FImpl>
std::vector<std::vector<Complex>> TDMixingTopD<FImpl>::contract_D(const typename TDMixingTopD<FImpl>::SlicedPropagator &half_if, const typename TDMixingTopD<FImpl>::SlicedPropagator &half_fi)
{
    Gamma G5(Gamma::Algebra::Gamma5);
    Gamma GT(Gamma::Algebra::GammaT);

    // Kept general in case anyone ever wants to play with this
    Gamma Gsrc = G5;
    Gamma Gsnk = Gsrc; // no conj on final interpolator for D-Dbar mixing

    int Nt = half_if.size();

    std::vector<std::vector<Complex>> corr(Nt, std::vector<Complex>(Nt, 0.));
    for (int t1 = 0; t1 < Nt; t1++)
    {
        for (int t2 = 0; t2 < Nt; t2++)
        {
            corr[t1][t2] = TensorRemove(trace(half_if[t1] * Gsrc * half_fi[t2] * Gsnk));
        }
    }

    return corr;
};

template <typename FImpl>
std::pair<typename TDMixingTopD<FImpl>::PropagatorField, typename TDMixingTopD<FImpl>::PropagatorField> TDMixingTopD<FImpl>::GH_VVAA_cap(const TDMixingTopD<FImpl>::PropagatorField &prop)
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

    SitePropagator spId = Zero();
    for (int s = 0; s < 4; s++)
    {
        for (int c = 0; c < 3; c++)
        {
            spId()(s, s)(c, c) = 1.;
        }
    }

    PropagatorField GTrPropG_VVAA(grid);
    GTrPropG_VVAA = Zero();
    PropagatorField GPropG_VVAA(grid);
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
void TDMixingTopD<FImpl>::setup(void)
{

    GridCartesian *grid = envGetGrid(FermionField);
    envTmp(std::vector<PropagatorField>, "GdsG_pp", 1, 2, PropagatorField(env().getGrid()));
    auto &ql1 = envGet(std::vector<PropagatorField *>, par().qLoop1);
    int Neta = ql1.size();
    const int Nt{env().getDim(Tdir)};
    envTmp(std::vector<SlicedPropagator>, "half_lr", 1, 2 * Neta, SlicedPropagator(Nt));
    envTmp(std::vector<SlicedPropagator>, "half_rl", 1, 2 * Neta, SlicedPropagator(Nt));
    envTmp(std::vector<std::vector<std::vector<Complex>>>, "buf", 1, Neta * Neta, std::vector<std::vector<Complex>>(Nt, std::vector<Complex>(Nt)));

    if (par().qLoop1 != par().qLoop2)
    {
        HADRONS_ERROR(Argument, "Current implementation for identical loops only");
    }
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology D" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    LOG(Message) << "qLoop1  : " << par().qLoop1 << std::endl;
    LOG(Message) << "qLoop2  : " << par().qLoop2 << std::endl;

    std::vector<Result> result;
    Result res;

    const int Nt{env().getDim(Tdir)};
    GridCartesian *grid = envGetGrid(FermionField);

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);
    auto &ql1 = envGet(std::vector<PropagatorField *>, par().qLoop1);
    auto &ql2 = envGet(std::vector<PropagatorField *>, par().qLoop2);

    // this is assuming both loops are identical
    int Neta = ql1.size();
    bool same_loop_noise = true;

    // std::vector<SlicedPropagator> half_lr(2 * Neta, SlicedPropagator(Nt));
    // std::vector<SlicedPropagator> half_rl(2 * Neta, SlicedPropagator(Nt));
    envGetTmp(std::vector<SlicedPropagator>, half_lr);
    envGetTmp(std::vector<SlicedPropagator>, half_rl);

    // parity +, parity -
    std::vector<Gamma> parityG = {Gamma(Gamma::Algebra::Identity), Gamma(Gamma::Algebra::Gamma5)};
    envGetTmp(std::vector<PropagatorField>, GdsG_pp);
    // std::vector<std::vector<std::vector<Complex>>> buf(Neta * Neta, std::vector<std::vector<Complex>>(Nt, std::vector<Complex>(Nt)));
    envGetTmp(std::vector<std::vector<std::vector<Complex>>>, buf);
    std::vector<std::vector<Complex>> tmp = std::vector<std::vector<Complex>>(Nt, std::vector<Complex>(Nt, 0.));
    for (int p = 0; p < 2; p++)
    {
        for (int i = 0; i < Neta; i++)
        {
            // here one has to add ql2 if one wants them to be allowed to be different
            auto tmp = GH_VVAA_cap(*ql1[i]);
            GdsG_pp[0] = tmp.first;  // r1
            GdsG_pp[1] = tmp.second; // r2
            for (int r = 0; r < 2; r++)
            {
                half_lr[i + Neta * r] = contract_D_half(qcl, qur, GdsG_pp[r] * parityG[p]);
                half_rl[i + Neta * r] = contract_D_half(qcr, qul, GdsG_pp[r] * parityG[p]);
            }
        }

        for (int r = 0; r < 2; r++)
        {
            for (int s = 0; s < 2; s++)
            {
                res.info.rr = std::to_string(r + 1) + std::to_string(s + 1);
                res.info.parity = (p == 0) ? "+" : "-";
                for (int i = 0; i < Neta; i++)
                {
                    for (int j = 0; j < Neta; j++)
                    {
                        buf[j + Neta * i] = contract_D(half_lr[i + Neta * r], half_rl[j + Neta * s]);
                    }
                }

                // Average noises up to imax (+ remove diagonal terms)
                for (int imax = 1; imax <= Neta; imax++)
                {
                    for (int t0 = 0; t0 < Nt; t0++)
                    {
                        std::fill(tmp[t0].begin(), tmp[t0].end(), 0.0);
                    }
                    const double norm = (!same_loop_noise) ? 1.0 / (imax * imax) : (imax > 1 ? 1.0 / (imax * (imax - 1)) : 1.0);
                    for (int i = 0; i < imax; i++)
                    {
                        for (int j = 0; j < imax; j++)
                        {
                            if (i == j and same_loop_noise)
                                continue;
                            const auto &c = buf[j + Neta * i];
                            for (int t1 = 0; t1 < Nt; t1++)
                            {
                                for (int t2 = 0; t2 < Nt; t2++)
                                {
                                    tmp[t1][t2] += c[t1][t2] * norm;
                                }
                            }
                        }
                    }
                    res.info.eta_max = std::to_string(imax);
                    res.corr.clear();
                    res.corr = tmp;
                    result.push_back(res);
                }
            }
        }
    }

    // save result, and hand it to environment
    saveResult(par().output, "DMixingTopD", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopD_hpp_
