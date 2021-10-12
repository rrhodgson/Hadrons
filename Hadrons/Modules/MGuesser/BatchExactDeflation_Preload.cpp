#include <Hadrons/Modules/MGuesser/BatchExactDeflation_Preload.hpp>

using namespace Grid;
using namespace Hadrons;
using namespace MGuesser;

template class Grid::Hadrons::MGuesser::TBatchExactDeflation_Preload<BaseFermionEigenPack<FIMPL>>;
template class Grid::Hadrons::MGuesser::TBatchExactDeflation_Preload<BaseFermionEigenPack<FIMPLF>>;
