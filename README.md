# boltzbub

…simple Lattice Boltzmann code experiments written in C++.

![channel flow visualized using ParaView](screenshot/channel_flow.png)

## Experiments

| Code                                  | Description                                                     |
| -                                     | -                                                               |
| `lid_driven_cavity.cc`                | Lid driven cavity using (moving wall) bounce back conditions    |
| `poiseuille.cc`                       | Poiseuille flow with velocity BC inflow and density BC outflow  |
| `channel.cc`                          | Channel flow with some obstacles and Dirichlet inflow condition |

## Build

	git clone https://github.com/KnairdA/boltzbub.git
	cd boltzbub
	nix-shell
	mkdir build
	cd build
	cmake ..
	make
