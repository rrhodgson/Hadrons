/*
 * DMixingTopA.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_DMixingTopA_hpp_
#define Hadrons_MContraction_DMixingTopA_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopA                                        *
 *                (Fig. 4 (A) in arxiv:2504.16189)
 *                   qCL        qInt1        qUR
 *                 /-->--\   /---->----\   /-->--\
 *                /       \ /           \ /       \
 *               /       ┌───┐         ┌───┐       \
 *           g5 *        | r |         | r'|        * g5
 *               \       └───┘         └───┘       /
 *                \       / \           / \       /
 *                 \--<--/   \----<----/   \--<--/
 *                   qUL        qInt2        qCR
 *         tsrc           t1            t2          tsnk
 *
 * Four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 *
 * p = +: GA x GB =   V x V + A x A
 * p = -: GA x GB = - A x V - V x A
 *
 * Contractions: [...] = tr(...)
 * rr' = 11:
 *  [qInt1 * GA1 * qCL * g5 * qUL * GB1 * qInt2 * GA2 * qCR * g5 * qUR * GB2]
 * rr' = 12:
 *  [qInt1 * GA1 * qCL * g5 * qUL * GB1 * qInt2 * GA2]*[qCR * g5 * qUR * GB2]
 * rr' = 21:
 *  [qCL * g5 * qUL * GB1]*[qInt2 * GA2 * qCR * g5 * qUR * GB2 * qInt1 * GA1]
 * rr' = 22:
 *  [qCL * g5 * qUL * GB1]*[qInt1 * GA1 * qInt2 * GA2]*[qCR * g5 * qUR * GB2]
 *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MContraction)

class DMixingTopAPar : Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(DMixingTopAPar,
                                    std::string, qULeft,
                                    std::string, qCLeft,
                                    std::string, qURight,
                                    std::string, qCRight,
                                    std::string, qInt,
                                    std::string, points,
                                    std::string, output);
};

template <typename FImpl>
class TDMixingTopA : public Module<DMixingTopAPar>
{
public:
    FERM_TYPE_ALIASES(FImpl, );
    class Metadata : Serializable
    {
    public:
        GRID_SERIALIZABLE_CLASS_MEMBERS(Metadata,
                                        int        , r,
                                        int        , s,
                                        std::string, parity);
    };
    typedef Correlator<Metadata, std::vector<Complex>> Result;

    using Parity   = typename DMixingUtils<FImpl>::Parity;
    using OpStruct = typename DMixingUtils<FImpl>::OpStruct;

    const std::array<Gamma,8> &GHs     = DMixingUtils<FImpl>::GHs;
    const std::array<Gamma,2> &parityG = DMixingUtils<FImpl>::parityG;

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
    virtual std::vector<std::vector<Complex>> contract_A(const PropagatorField &GcuGl, const PropagatorField &GcuGr, const std::vector<Coordinate *> &xs, const std::vector<PropagatorField *> &ds_prop_pt);
};

MODULE_REGISTER_TMP(DMixingTopA, TDMixingTopA<FIMPL>, MContraction);

/******************************************************************************
 *                        TDMixingTopA implementation                         *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TDMixingTopA<FImpl>::TDMixingTopA(const std::string name)
    : Module<DMixingTopAPar>(name)
{
}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TDMixingTopA<FImpl>::getInput(void)
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
std::vector<std::vector<Complex>> TDMixingTopA<FImpl>::contract_A(
    const PropagatorField &GcuG_l,
    const PropagatorField &GcuG_r,
    const std::vector<Coordinate *> &xs,
    const std::vector<PropagatorField *> &ds_prop_pt)
{
    int Nt = env().getDim(Tdir);
    Gamma g5(Gamma::Algebra::Gamma5);
    std::vector<std::vector<Complex>> corr(Nt, std::vector<Complex>(Nt));

    SlicedPropagator B;
    B.reserve(Nt);

    for (int t1 = 0; t1 < Nt; t1++)
    {
        startTimer("peekSite");
        const auto A = peekSite(GcuG_l, *xs[t1]);
        stopTimer("peekSite");

        startTimer("mult");
        const PropagatorField &ds = *ds_prop_pt[t1];
        const PropagatorField dsD = g5 * adj(ds) * g5;
        PropagatorField tmp = dsD * GcuG_r * ds;
        stopTimer("mult");

        startTimer("sliceSum");
        sliceSum(tmp, B, Tp);
        stopTimer("sliceSum");

        startTimer("trace");
        for (int t2 = 0; t2 < Nt; t2++)
            corr[t1][t2] = TensorRemove(trace(A * B[t2]));
        stopTimer("trace");
    }

    return corr;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopA<FImpl>::setup(void)
{
    envTmpLat(PropagatorField, "qcul");
    envTmpLat(PropagatorField, "qcur");
    envTmp(std::vector<PropagatorField>, "GcuG_l", 1, 2, PropagatorField(env().getGrid()));
    envTmp(std::vector<PropagatorField>, "GcuG_r", 1, 2, PropagatorField(env().getGrid()));
    envTmpLat(PropagatorField, "half_l");
    envTmpLat(PropagatorField, "half_r");

    envCreate(HadronsSerializable, getName(), 1, 0);
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopA<FImpl>::execute(void)
{
    LOG(Message) << "Computing D-meson mixing diagram, topology A" << std::endl;
    LOG(Message) << "qULeft  : " << par().qULeft << std::endl;
    LOG(Message) << "qCLeft  : " << par().qCLeft << std::endl;
    LOG(Message) << "qURight : " << par().qURight << std::endl;
    LOG(Message) << "qCRight : " << par().qCRight << std::endl;
    LOG(Message) << "qInt    : " << par().qInt << std::endl;
    LOG(Message) << "points  : " << par().points << std::endl;

    std::vector<Result> result;

    const int Nt{env().getDim(Tdir)};

    auto &qul = envGet(PropagatorField, par().qULeft);
    auto &qcl = envGet(PropagatorField, par().qCLeft);
    auto &qur = envGet(PropagatorField, par().qURight);
    auto &qcr = envGet(PropagatorField, par().qCRight);
    auto &qin = envGet(std::vector<PropagatorField *>, par().qInt);
    auto &pts = envGet(std::vector<Coordinate *>, par().points);

    envGetTmp(PropagatorField, qcul);
    envGetTmp(PropagatorField, qcur);
    envGetTmp(std::vector<PropagatorField>, GcuG_l);
    envGetTmp(std::vector<PropagatorField>, GcuG_r);
    envGetTmp(PropagatorField, half_l);
    envGetTmp(PropagatorField, half_r);

    Gamma g5(Gamma::Algebra::Gamma5);

    qcul = qcl * adj(qul) * g5; // qcl * g5 * (g5 * adj(qul) * g5)
    qcur = qcr * adj(qur) * g5; // qcr * g5 * (g5 * adj(qur) * g5)

    startTimer("GH_cap");
    DMixingUtils<FImpl>::GH_cap(GcuG_l, qcul, Parity::Pos);
    DMixingUtils<FImpl>::GH_cap(GcuG_r, qcur, Parity::Pos);
    stopTimer("GH_cap");

    for (const auto p : {Parity::Pos,Parity::Neg})
    {
        for (const auto r : {OpStruct::One,OpStruct::Two})
        {
            half_l = parityG[p] * GcuG_l[r];

            for (const auto s : {OpStruct::One,OpStruct::Two})
            {
                half_r = parityG[p] * GcuG_r[s];

                Result res;
                res.info.parity = DMixingUtils<FImpl>::toString(p);
                res.info.r      = DMixingUtils<FImpl>::toInt(r);
                res.info.s      = DMixingUtils<FImpl>::toInt(s);

                res.corr = contract_A(half_l, half_r, pts, qin);

                result.push_back(res);
            }
        }
    }

    // save result, and hand it to environment
    saveResult(par().output, "DMixingTopA", result);
    auto &out = envGet(HadronsSerializable, getName());
    out = result;
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MContraction_DMixingTopA_hpp_
