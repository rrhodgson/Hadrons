#include <Hadrons/Modules/MGuesser/CoarseExactDeflation.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPL,HADRONS_DEFAULT_LANCZOS_NBASIS>>;
template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPLF,HADRONS_DEFAULT_LANCZOS_NBASIS>>;

template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPL,30>>;
template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPLF,30>>;

template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPL,150>>;
template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPLF,150>>;

template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPL,400>>;
template class Grid::Hadrons::MGuesser::TCoarseExactDeflation<CoarseFermionEigenPack<FIMPLF,400>>;

