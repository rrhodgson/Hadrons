/*
 * DMixingTopD.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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
#ifndef Hadrons_MContraction_DMixingTopD_hpp_
#define Hadrons_MContraction_DMixingTopD_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>
#include <Hadrons/Serialization.hpp>
#include <Hadrons/Modules/MContraction/DMixingUtils.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                         DMixingTopD                                        *
 *                    (Fig. 4 (D) in arxiv:2504.16189)
 *             qCL   ┌───┐                 qUR
 *         /---->----| r |----------------->--------------\
 *        /          └───┘                                 \
 *       /          /     \          qLoop2 /-<-\           \
 *   g5 *           \     /                /     \           * g5
 *       \           \->-/ qLoop1          \     /          /
 *        \                                 ┌───┐          /
 *         \--------------<-----------------| r'|----<----/
 *                       qUL                └───┘  qCR
 *  tsrc               t1                    t2              tsnk
 *
 * Four configurations for the two weak Hamiltonians M_r M_{r'}
 * (cf. Fig. 3 in arxiv:2504.16189)
 *
 * p = +: GA x GB =   V x V + A x A
 * p = -: GA x GB = - A x V - V x A
 *
 * Contractions: [...] = tr(...)
 * rr' = 11:
 *  [qCR * g5 * qUR * GB1 * qLoop1 * GA1 * qCL * g5 * gUL * GB2 * gLoop2 * GA2]
 * rr' = 12:
 *  [qCR * g5 * qUR * GB1 * qLoop1 * GA1 * qCL * g5 * gUL * GB2]*[gLoop2 * GA2]
 * rr' = 21:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * gUL * GB2 * gLoop2 * GA2]*[qLoop1 * GA1]
 * rr' = 22:
 *  [qCR * g5 * qUR * GB1 * qCL * g5 * gUL * GB2]*[gLoop2 * GA2]*[qLoop1 * GA1]
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
    virtual void contract_D_half(const PropagatorField &prop_c,
                                 const PropagatorField &prop_u,
                                 const PropagatorField &loop,
                                 SlicedPropagator &out);
    virtual void contract_D(const SlicedPropagator &half_lr,
                           const SlicedPropagator &half_rl,
                           std::vector<std::vector<Complex>> &sum);
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
inline void TDMixingTopD<FImpl>::contract_D_half(
    const typename TDMixingTopD<FImpl>::PropagatorField &prop_c,
    const typename TDMixingTopD<FImpl>::PropagatorField &prop_u,
    const typename TDMixingTopD<FImpl>::PropagatorField &loop,
    SlicedPropagator &out)
{
    Gamma g5(Gamma::Algebra::Gamma5);

    PropagatorField tmp = g5 * adj(prop_u) * g5 * loop * prop_c;
    sliceSum(tmp, out, Tp);
};

template <typename FImpl>
void TDMixingTopD<FImpl>::contract_D(
    const typename TDMixingTopD<FImpl>::SlicedPropagator &half_lr,
    const typename TDMixingTopD<FImpl>::SlicedPropagator &half_rl,
    std::vector<std::vector<Complex>> &sum)
{
    Gamma g5(Gamma::Algebra::Gamma5);
    const int Nt = (int)half_lr.size();

    for (int t1 = 0; t1 < Nt; t1++)
        for (int t2 = 0; t2 < Nt; t2++)
            sum[t1][t2] += TensorRemove(trace(half_lr[t1] * g5 * half_rl[t2] * g5));
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TDMixingTopD<FImpl>::setup(void)
{

    GridCartesian *grid = envGetGrid(FermionField);

    envTmpLat(PropagatorField, "half_lr_i");
    envTmpLat(PropagatorField, "half_rl_i");
    envTmpLat(PropagatorField, "half_lr_j");
    envTmpLat(PropagatorField, "half_rl_j");
    envTmp(std::vector<PropagatorField>, "GdsG", 1, 2, PropagatorField(env().getGrid()));

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

    envGetTmp(SlicedPropagator, half_lr_i);
    envGetTmp(SlicedPropagator, half_rl_i);
    envGetTmp(SlicedPropagator, half_lr_j);
    envGetTmp(SlicedPropagator, half_rl_j);
    envGetTmp(std::vector<PropagatorField>, GdsG);

    std::vector<std::vector<Complex>> tmpRes = std::vector<std::vector<Complex>>(Nt, std::vector<Complex>(Nt, 0.));
    std::vector<std::vector<Complex>> tmpSum = std::vector<std::vector<Complex>>(Nt, std::vector<Complex>(Nt, 0.));

    for (int p = 0; p < 2; p++)
    {
        for (int r = 0; r < 2; r++)
        {
            for (int s = 0; s < 2; s++)
            {
                for (auto &row : tmpSum)
                    std::fill(row.begin(), row.end(), Complex(0.0));

                for (int i = 0; i < Neta; i++) // imax = i + 1
                {
                    DMixingUtils<FImpl>::GH_cap(GdsG, *ql1[i], p);

                    contract_D_half(qcl, qur, GdsG[r], half_lr_i);
                    contract_D_half(qcr, qul, GdsG[s], half_rl_i);

                    for (int j = 0; j < i; j++)
                    {
                        DMixingUtils<FImpl>::GH_cap(GdsG, *ql1[j], p);
                        contract_D_half(qcl, qur, GdsG[r], half_lr_j);
                        contract_D_half(qcr, qul, GdsG[s], half_rl_j);

                        contract_D(half_lr_i, half_rl_j, tmpSum);
                        contract_D(half_lr_j, half_rl_i, tmpSum);
                    }
                    if (!same_loop_noise)
                    {
                        // include diagonal term when loops use different noises
                        contract_D(half_lr_i, half_rl_i, tmpSum);
                    }

                    const double norm = (!same_loop_noise)
                                            ? double((i + 1) * (i + 1))
                                            : (i > 0 ? double((i + 1) * i) : 1.0);

                    for (int t1 = 0; t1 < Nt; t1++)
                        for (int t2 = 0; t2 < Nt; t2++)
                            tmpRes[t1][t2] = tmpSum[t1][t2] / norm;

                    Result res;
                    res.info.parity = (p == 0) ? "+" : "-";
                    res.info.rr = std::to_string(r + 1) + std::to_string(s + 1);
                    res.info.eta_max = std::to_string(i + 1);
                    res.corr = tmpRes;
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
