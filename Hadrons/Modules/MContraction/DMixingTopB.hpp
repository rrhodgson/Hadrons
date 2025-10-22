/*
 * DMixingTopB.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2025
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Matteo Di Carlo <matteo.dicarlo@cern.ch>
 * Author: Felix Erben <felix.erben@cern.ch>
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
 *                 /---->----| r |--------->-----------\
 *                /          └───┘--<--\                \
 *               /              \        \ qInt2         \
 *           g5 *                \        \               * g5
 *               \          qInt1 \        \             /
 *                \                \-->--┌───┐          /
 *                 \---------<-----------| r'|----<----/
 *                          qUL          └───┘  qCR
 *          tsrc               t1          t2             tsnk
 *
 * Four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 *
 * p = +: GA x GB =   V x V + A x A
 * p = -: GA x GB = - A x V - V x A
 *
 * Contractions: [...] = tr(...)
 * rr' = 11:
 *  [qInt1 * GA1 * qCL * g5 * qUL * GB2]*[qCR * g5 * qUR * GB1 * qInt2 * GA2]
 * rr' = 12:
 *  [qInt1 * GA1 * qCL * g5 * qUL * GB2 * qCR * g5 * qUR * GB1 * qInt2 * GA2]
 * rr' = 21:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * qUL * GB2 * qInt1 * GA1 * qInt2 * GA2]
 * rr' = 22:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * qUL * GB2]*[qInt1 * GA1 * qInt2 * GA2]
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

    using Parity   = typename DMixingUtils<FImpl>::Parity;
    using OpStruct = typename DMixingUtils<FImpl>::OpStruct;

    const std::array<Gamma,8> &GHs     = DMixingUtils<FImpl>::GHs;
    const std::array<Gamma,2> &parityG = DMixingUtils<FImpl>::parityG;

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
    virtual std::vector<std::vector<Complex>> contract_B(
        const PropagatorField &ci,
        const PropagatorField &ui,
        const PropagatorField &cf,
        const PropagatorField &uf,
        const std::vector<PropagatorField *> &ds_prop_pt,
        const std::vector<Coordinate *> &xs,
        const Parity p,
        const OpStruct rL,
        const OpStruct rR);
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
std::vector<std::vector<Complex>> TDMixingTopB<FImpl>::contract_B(
    const PropagatorField &ci,
    const PropagatorField &ui,
    const PropagatorField &cf,
    const PropagatorField &uf,
    const std::vector<PropagatorField *> &ds_prop_pt,
    const std::vector<Coordinate *> &xs,
    const Parity p,
    const OpStruct rL,
    const OpStruct rR)
{
    int Nt = env().getDim(Tdir);
    Gamma g5(Gamma::Algebra::Gamma5);

    std::vector<std::vector<Complex>> corr(Nt, std::vector<Complex>(Nt, 0.));
    std::vector<TComplex> buf;
    buf.reserve(Nt);

    const Gamma parity = parityG[p];

    const bool same_r = (rR == rL);

    for (int t1 = 0; t1 < Nt; t1++)
    {
        const PropagatorField &ds = *ds_prop_pt[t1];
        const PropagatorField dsD = g5 * adj(ds) * g5;

        const PropagatorField cui = peekSite(ci, *xs[t1]) * adj(ui) * g5;
        const PropagatorField cuf = cf * adj(peekSite(uf, *xs[t1])) * g5;

        // product of traces
        auto tr_same = [&](const auto &trA, const auto &trB)
        {
            for (const auto &GH2 : GHs)
            {
                startTimer("trace");
                LatticeComplex tmp = trace(trA * GH2) * trace(trB * parity * GH2);
                stopTimer("trace");

                startTimer("sliceSum");
                sliceSum(tmp, buf, Tp);
                stopTimer("sliceSum");

                auto &row = corr[t1];
                for (int t2 = 0; t2 < Nt; t2++)
                    row[t2] += TensorRemove(buf[t2]);
            }
        };
        // trace of product
        auto tr_diff = [&](const auto &trA, const auto &trB)
        {
            for (const auto &GH2 : GHs)
            {
                startTimer("trace");
                LatticeComplex tmp = trace(trA * GH2 * trB * parity * GH2);
                stopTimer("trace");

                startTimer("sliceSum");
                sliceSum(tmp, buf, Tp);
                stopTimer("sliceSum");

                auto &row = corr[t1];
                for (int t2 = 0; t2 < Nt; t2++)
                    row[t2] += TensorRemove(buf[t2]);
            }
        };

        switch (rL)
        {
            case OpStruct::One:
                for (const auto &GH1 : GHs)
                {
                    startTimer("trA & trB");
                    PropagatorField trA = ds * parity * GH1 * cui;
                    PropagatorField trB = cuf * GH1 * dsD;
                    stopTimer("trA & trB");
                    same_r ? tr_same(trA, trB) : tr_diff(trA, trB);
                }
                break;

            case OpStruct::Two:
                for (const auto &GH1 : GHs)
                {
                    startTimer("trA & trB");
                    PropagatorField trA = cuf * GH1 * cui;
                    PropagatorField trB = ds * parity * GH1 * dsD;
                    stopTimer("trA & trB");
                    same_r ? tr_same(trA, trB) : tr_diff(trA, trB);
                }
                break;

            default:
                HADRONS_ERROR(Argument, "DMixingTopB: Invalid RS value");
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

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);
    auto &qin = envGet(std::vector<PropagatorField *>, par().qInt);
    auto &pts = envGet(std::vector<Coordinate *>, par().points);

    for (const auto p : {Parity::Pos,Parity::Neg})
    {
        for (const auto r : {OpStruct::One,OpStruct::Two})
        {
            for (const auto s : {OpStruct::One,OpStruct::Two})
            {
                Result res;
                res.info.parity = (p == Parity::Pos) ? "+" : "-";
                res.info.rr = std::to_string(r + 1) + std::to_string(s + 1);

                res.corr = contract_B(qcl, qul, qcr, qur, qin, pts, p, r, s);
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
