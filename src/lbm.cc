#include "lbm.h"

void Cell::equilibrize(Density d, Velocity v) {
	for ( int i = -1; i <= 1; ++i ) {
		for ( int j = -1; j <= 1; ++j ) {
			get(i,j) = weight.get(i,j) * d * (1 + 3*v.comp(i,j) + 4.5*sq(v.comp(i,j)) - 1.5*sq(v.norm()));
		}
	}
}

double Cell::sum() const {
	return get(-1, 1) + get( 0, 1) + get( 1, 1) + get(-1, 0) + get( 0, 0) + get( 1, 0) + get(-1,-1) + get( 0,-1) + get( 1,-1);
}

Velocity Cell::velocity(Density d) const {
	return 1./d * Velocity{
		get( 1, 0) - get(-1, 0) + get( 1, 1) - get(-1,-1) + get( 1,-1) - get(-1,1),
		get( 0, 1) - get( 0,-1) + get( 1, 1) - get(-1,-1) - get( 1,-1) + get(-1,1)
	};
}

void streamFluidCell(DataCellBuffer& pop, std::size_t x, std::size_t y) {
	if ( x != 0 && x != pop.dimX()-1 && y != 0 && y != pop.dimY()-1 ) {
		// stream internal cells
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				pop.curr(x+i,y+j).get(i,j) = pop.prev(x,y).get(i,j);
			}
		}
	} else {
		// stream boundary cells,
		// missing populations to be determined by boundary conditions
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				if ( (x > 0 || i >= 0) && x+i <= pop.dimX()-1 && (y > 0 || j >= 0) && y+j <= pop.dimY()-1 ) {
					pop.curr(x+i,y+j).get(i,j) = pop.prev(x,y).get(i,j);
				}
			}
		}
	}
}

void collideFluidCell(
	double omega,
	DataCellBuffer& pop, FluidBuffer& fluid,
	std::size_t x, std::size_t y) {
	// compute equilibrium
	Cell eq;
	eq.equilibrize(fluid.density(x,y), fluid.velocity(x,y));

	// collide (BGK, relax towards equilibrium)
	if ( x != 0 && x != pop.dimX()-1 && y != 0 && y != pop.dimY()-1 ) {
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				pop.curr(x,y).get(i,j) = pop.curr(x,y).get(i,j) + omega * (eq.get(i,j) - pop.curr(x,y).get(i,j));
			}
		}
	} else {
		// partial collide for boundary cells
		for ( int i = -1; i <= 1; ++i ) {
			for ( int j = -1; j <= 1; ++j ) {
				if ( (x > 0 || i >= 0) && x+i <= pop.dimX()-1 && (y > 0 || j >= 0) && y+j <= pop.dimY()-1 ) {
					pop.curr(x,y).get(i,j) = pop.curr(x,y).get(i,j) + omega * (eq.get(i,j) - pop.curr(x,y).get(i,j));
				}
			}
		}
	}
}
