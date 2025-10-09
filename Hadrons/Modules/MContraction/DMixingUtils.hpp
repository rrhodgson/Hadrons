/*
 * DMixingUtils.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
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

#ifndef Hadrons_MContraction_DMixingUtils_hpp_
#define Hadrons_MContraction_DMixingUtils_hpp_

#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE
BEGIN_MODULE_NAMESPACE(MContraction)

template <typename FImpl>
class DMixingUtils
{
public:
    FERM_TYPE_ALIASES(FImpl, )
    static void GH_cap(std::vector<PropagatorField> &out, const PropagatorField &prop, const Gamma &parityG);
};

template <typename FImpl>
void DMixingUtils<FImpl>::GH_cap(
    std::vector<typename DMixingUtils<FImpl>::PropagatorField> &out,
    const typename DMixingUtils<FImpl>::PropagatorField &prop,
    const Gamma &parityG)
{
    assert(out.size() == 2);

    std::array<Gamma, 8> GHs{Gamma(Gamma::Algebra::GammaX),
                             Gamma(Gamma::Algebra::GammaY),
                             Gamma(Gamma::Algebra::GammaZ),
                             Gamma(Gamma::Algebra::GammaT),
                             Gamma(Gamma::Algebra::GammaXGamma5),
                             Gamma(Gamma::Algebra::GammaYGamma5),
                             Gamma(Gamma::Algebra::GammaZGamma5),
                             Gamma(Gamma::Algebra::GammaTGamma5)};

    SitePropagator spId(1.0);
    out[0] = Zero();
    out[1] = Zero();

    for (const auto &GH : GHs)
    {
        out[0] += GH * prop * parityG * GH;
        out[1] += spId * GH * trace(prop * parityG * GH);
    }
};

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif
