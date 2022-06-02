/*
 * LoadCoarseEigenPack.cpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
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
#include <Hadrons/Modules/MIO/LoadCoarseEigenPack.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MIO;

template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>,GIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF,HADRONS_DEFAULT_LANCZOS_NBASIS>,GIMPLF>;
template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS, FIMPLF>,GIMPL>;
#endif

template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL,150>,GIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF,150>,GIMPLF>;
template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL,150, FIMPLF>,GIMPL>;
#endif

template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL,400>,GIMPL>;
#ifdef GRID_DEFAULT_PRECISION_DOUBLE
template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPLF,400>,GIMPLF>;
template class Grid::Hadrons::MIO::TLoadCoarseEigenPack<CoarseFermionEigenPack<FIMPL,400, FIMPLF>,GIMPL>;
#endif
