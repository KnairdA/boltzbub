#include "boundary_conditions.h"

std::pair<Vector<int>, Vector<int>> neighbors(Vector<int> v) {
	if ( v[0] == 0 ) {
		return {
			{ -1, v[1] },
			{  1, v[1] }
		};
	} else if ( v[1] == 0 ) {
		return {
			{ v[0], -1 },
			{ v[0],  1 }
		};
	} else {
		return {
			{    0, v[1] },
			{ v[0], 0    }
		};
	}
}

void computeWallCell(DataCellBuffer& pop, Vector<std::size_t> cell, Vector<int> normal) {
	const auto [neighborA, neighborB] = neighbors(normal);

	pop.curr(cell).get(neighborA) = pop.curr(cell).get(-neighborA);
	pop.curr(cell).get(normal   ) = pop.curr(cell).get(-normal   );
	pop.curr(cell).get(neighborB) = pop.curr(cell).get(-neighborB);
}

void computeZouHeVelocityWallCell(DataCellBuffer& pop, Vector<std::size_t> cell, Vector<int> normal, double vX) {
	const auto [neighborA, neighborB] = neighbors(normal);

	const double rho = pop.curr(cell).get(-1,0) + pop.curr(cell).get(0,0) + pop.curr(cell).get(1,0)
		+ 2.*(
			pop.curr(cell).get(-neighborA) +
			pop.curr(cell).get(-normal   ) +
			pop.curr(cell).get(-neighborB)
		);

	pop.curr(cell).get(neighborA) = pop.curr(cell).get(-neighborA)
	                              + 0.5*( pop.curr(cell).get( 1,0) - pop.curr(cell).get(-1,0) - vX*rho );
	pop.curr(cell).get(normal   ) = pop.curr(cell).get(-normal   );
	pop.curr(cell).get(neighborB) = pop.curr(cell).get(-neighborB)
	                              + 0.5*( pop.curr(cell).get(-1,0) - pop.curr(cell).get( 1,0) + vX*rho );
}