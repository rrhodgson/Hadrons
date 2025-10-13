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

class DMixingTopBPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopBPar,
                                    std::string, qULeft,
                                    std::string, qCLeft,
                                    std::string, qURight,
                                    std::string, qCRight,
                                    std::string, qInt,
                                    std::string, points,
                                    std::string, output);
};

template <typename FImpl>
class TDMixingTopB : public Module<DMixingTopBPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        std::string, rr,
                                        std::string, parity);
    };
    typedef Correlator<Metadata, std::vector<Complex>> Result;

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
    // bespoke subcontractions
    virtual std::vector<std::vector<Complex>> contractB(const PropagatorField &ci, const PropagatorField &ui, const PropagatorField &cf, const PropagatorField &uf, const std::vector<PropagatorField*> &ds_prop_pt, const std::vector<Coordinate*> &xs, const std::array<Gamma, 8> &GHs, const Gamma &parityG, int rL, int rR);
};

MODULE_REGISTER_TMP(DMixingTopB, TDMixingTopB<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopB implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopB<FImpl>::TDMixingTopB(const std::string name)
    : Module<DMixingTopBPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopB<FImpl>::getInput(void)
{
    std::vector<std::string> in = {par().qULeft,
                                   par().qCLeft,
                                   par().qURight,
                                   par().qCRight,
                                   par().qInt,
                                   par().points};
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

template <typename FImpl>
std::vector<std::vector<Complex>> TDMixingTopB<FImpl>::contractB(
    const TDMixingTopB<FImpl>::PropagatorField &ci,
    const TDMixingTopB<FImpl>::PropagatorField &ui,
    const TDMixingTopB<FImpl>::PropagatorField &cf,
    const TDMixingTopB<FImpl>::PropagatorField &uf,
    const std::vector<typename TDMixingTopB<FImpl>::PropagatorField*> &ds_prop_pt,
    const std::vector<Coordinate*> &xs,
    const std::array<Gamma, 8> &GHs,
    const Gamma &parityG,
    const int rL,
    const int rR)
{
    int Nt = ci.Grid()->_fdimensions[3];
    Gamma g5(Gamma::Algebra::Gamma5);

    std::vector<std::vector<Complex>> corr(Nt, std::vector<Complex>(Nt, 0.));
    std::vector<TComplex> buf;

    const bool same_r = rR == rL;

    for (int t1 = 0; t1 < Nt; t1++)
    {

        const auto &ds = *ds_prop_pt[t1];
        const PropagatorField dsD = g5 * adj(ds) * g5;

        const PropagatorField cui = peekSite(ci, *xs[t1]) * adj(ui) * g5;
        const PropagatorField cuf = cf * adj(peekSite(uf, *xs[t1])) * g5;

        for (const auto &GH1 : GHs)
        {
            const PropagatorField trA = (rL == 0)
                                            ? PropagatorField(ds * parityG * GH1 * cui)
                                            : PropagatorField(cuf * GH1 * cui);

            const PropagatorField trB = (rL == 0)
                                            ? PropagatorField(cuf * GH1 * dsD)
                                            : PropagatorField(ds * parityG * GH1 * dsD);

            if (same_r) // product of traces
            {
                for (const auto &GH2 : GHs)
                {
                    LatticeComplex tmp = trace(trA * GH2) * trace(trB * parityG * GH2);
                    sliceSum(tmp, buf, Tp);
                    for (int t2 = 0; t2 < Nt; ++t2)
                        corr[t1][t2] += TensorRemove(buf[t2]);
                }
            }
            else // trace of product
            {
                for (const auto &GH2 : GHs)
                {
                    LatticeComplex tmp = trace(trA * GH2 * trB * parityG * GH2);
                    sliceSum(tmp, buf, Tp);
                    for (int t2 = 0; t2 < Nt; ++t2)
                        corr[t1][t2] += TensorRemove(buf[t2]);
                }
            }
        }
    }

    return corr;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopB<FImpl>::setup(void)
{
    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopB<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology B" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    LOG(Message) << "qInt    : " << par().qInt << std::endl;
    LOG(Message) << "points  : " << par().points << std::endl;

    std::vector<Result> result;
    Result res;

    const int Nt{env().getDim(Tdir)};

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);
    auto &qi = envGet(std::vector<PropagatorField*>, par().qInt);
    auto &points = envGet(std::vector<Coordinate*>, par().points);

    std::array<Gamma, 8> GHs{Gamma(Gamma::Algebra::GammaX),
                             Gamma(Gamma::Algebra::GammaY),
                             Gamma(Gamma::Algebra::GammaZ),
                             Gamma(Gamma::Algebra::GammaT),
                             Gamma(Gamma::Algebra::GammaXGamma5),
                             Gamma(Gamma::Algebra::GammaYGamma5),
                             Gamma(Gamma::Algebra::GammaZGamma5),
                             Gamma(Gamma::Algebra::GammaTGamma5)};

    std::vector<Gamma> parityG =
        {Gamma(Gamma::Algebra::Identity), Gamma(Gamma::Algebra::Gamma5)};

    for (int p = 0; p < 2; p++)
    {
        for (int r = 0; r < 2; r++)
        {
            for (int s = 0; s < 2; s++)
            {
                res.info.rr = std::to_string(r + 1) + std::to_string(s + 1);
                res.info.parity = (p == 0) ? "+" : "-";
                res.corr.clear();
                res.corr =
                    contractB(qcl, qul, qcr, qur, qi, points, GHs, parityG[p], r, s);
                result.push_back(res);
            }
        }
    }

    // save result, and hand it to environment
    saveResult(par().output, "DMixingTopB", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopB_hpp_
