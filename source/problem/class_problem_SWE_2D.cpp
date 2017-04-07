#include "../general_definitions.h"
#include "class_problem_SWE_2D.h"

PROBLEM::PROBLEM() {
    this->mesh = new MESH(2,0);

    if(!(this->mesh->interfaces.empty())){
        if (this->mesh->interfaces.find(INTERNAL) != this->mesh->interfaces.end()) {
            this->internal_interfaces = this->mesh->interfaces.find(INTERNAL)->second;
        }
        if (this->mesh->interfaces.find(OCEAN) != this->mesh->interfaces.end()) {
            this->ocean_interfaces = this->mesh->interfaces.find(OCEAN)->second;
        }
        if (this->mesh->interfaces.find(LAND) != this->mesh->interfaces.end()) {
            this->land_interfaces = this->mesh->interfaces.find(LAND)->second;
        }
        if (this->mesh->interfaces.find(FLOW) != this->mesh->interfaces.end()) {
            this->flow_interfaces = this->mesh->interfaces.find(FLOW)->second;
        }
    }
    else {
        printf("\n");
        printf("PROBLEM - Fatal error!\n");
        printf("NO INTERFACES");
        exit(1);
    }
}

PROBLEM::~PROBLEM() {
    delete this->mesh;
}

void PROBLEM::EETimeStepper(int n_steps) {
    for (int step = 0; step < n_steps; step++) {
        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u[UA][i] += dt*it->second->u[D_UA][i];
                it->second->u[VA][i] += dt*it->second->u[D_VA][i];
                it->second->u[H][i] += dt*it->second->u[D_H][i];
            }
        }
    }
}

void PROBLEM::RK2TimeStepper(int n_steps) {
    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        it->second->u_substep.push_back(new double*[3]);
        for (int j = 0; j < 3; j++) {
            it->second->u_substep.at(0)[j] = new double[it->second->number_bf];
        }
    }

    for (int step = 0; step < n_steps; step++) {
        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u_substep.at(U0)[0][i] = it->second->u[UA][i];
                it->second->u_substep.at(U0)[1][i] = it->second->u[VA][i];
                it->second->u_substep.at(U0)[2][i] = it->second->u[H][i];

                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + (dt / 2.0)*it->second->u[D_UA][i];
                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + (dt / 2.0)*it->second->u[D_VA][i];
                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + (dt / 2.0)*it->second->u[D_H][i];
            }
        }

        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + dt*it->second->u[D_UA][i];
                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + dt*it->second->u[D_VA][i];
                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + dt*it->second->u[D_H][i];
            }
        }
    }

    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        for (int j = 0; j < 3; j++) {
            delete[] it->second->u_substep.at(0)[j];
        }
        delete[] it->second->u_substep.at(0);
        it->second->u_substep.clear();
    }
}

void PROBLEM::RK3TimeStepper(int n_steps) {
    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        for (int i = 0; i < 3; i++) {
            it->second->u_substep.push_back(new double*[3]);
            for (int j = 0; j < 3; j++) {
                it->second->u_substep.at(i)[j] = new double[it->second->number_bf];
            }
        }
    }

    for (int step = 0; step < n_steps; step++) {
        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u_substep.at(U0)[0][i] = it->second->u[UA][i];
                it->second->u_substep.at(U0)[1][i] = it->second->u[VA][i];
                it->second->u_substep.at(U0)[2][i] = it->second->u[H][i];

                it->second->u_substep.at(K1)[0][i] = it->second->u[D_UA][i];
                it->second->u_substep.at(K1)[1][i] = it->second->u[D_VA][i];
                it->second->u_substep.at(K1)[2][i] = it->second->u[D_H][i];

                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + (dt / 2.0)*it->second->u_substep.at(K1)[0][i];
                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + (dt / 2.0)*it->second->u_substep.at(K1)[1][i];
                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + (dt / 2.0)*it->second->u_substep.at(K1)[2][i];
            }
        }

        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u_substep.at(K2)[0][i] = it->second->u[D_UA][i];
                it->second->u_substep.at(K2)[1][i] = it->second->u[D_VA][i];
                it->second->u_substep.at(K2)[2][i] = it->second->u[D_H][i];

                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + dt*
                    (-it->second->u_substep.at(K1)[0][i] + 2 * it->second->u_substep.at(K2)[0][i]);
                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + dt*
                    (-it->second->u_substep.at(K1)[1][i] + 2 * it->second->u_substep.at(K2)[1][i]);
                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + dt*
                    (-it->second->u_substep.at(K1)[2][i] + 2 * it->second->u_substep.at(K2)[2][i]);
            }
        }

        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + (dt / 6.0)*
                    (it->second->u_substep.at(K1)[0][i] + 4 * it->second->u_substep.at(K2)[0][i] +
                        it->second->u[D_UA][i]);

                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + (dt / 6.0)*
                    (it->second->u_substep.at(K1)[1][i] + 4 * it->second->u_substep.at(K2)[1][i] +
                        it->second->u[D_VA][i]);

                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + (dt / 6.0)*
                    (it->second->u_substep.at(K1)[2][i] + 2 * it->second->u_substep.at(K2)[2][i] +
                        it->second->u[D_H][i]);
            }
        }
    }

    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                delete[] it->second->u_substep.at(i)[j];
            }
            delete[] it->second->u_substep.at(i);
        }
        it->second->u_substep.clear();
    }
}

void PROBLEM::RK4TimeStepper(int n_steps) {
    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        for (int i = 0; i < 4; i++) {
            it->second->u_substep.push_back(new double*[3]);
            for (int j = 0; j < 3; j++) {
                it->second->u_substep.at(i)[j] = new double[it->second->number_bf];
            }
        }
    }

    for (int step = 0; step < n_steps; step++) {
        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u_substep.at(U0)[0][i] = it->second->u[UA][i];
                it->second->u_substep.at(U0)[1][i] = it->second->u[VA][i];
                it->second->u_substep.at(U0)[2][i] = it->second->u[H][i];

                it->second->u_substep.at(K1)[0][i] = it->second->u[D_UA][i];
                it->second->u_substep.at(K1)[1][i] = it->second->u[D_VA][i];
                it->second->u_substep.at(K1)[2][i] = it->second->u[D_H][i];

                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + (dt / 2.0)*it->second->u_substep.at(K1)[0][i];
                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + (dt / 2.0)*it->second->u_substep.at(K1)[1][i];
                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + (dt / 2.0)*it->second->u_substep.at(K1)[2][i];
            }
        }

        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u_substep.at(K2)[0][i] = it->second->u[D_UA][i];
                it->second->u_substep.at(K2)[1][i] = it->second->u[D_VA][i];
                it->second->u_substep.at(K2)[2][i] = it->second->u[D_H][i];

                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + (dt / 2.0)*it->second->u_substep.at(K2)[0][i];
                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + (dt / 2.0)*it->second->u_substep.at(K2)[1][i];
                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + (dt / 2.0)*it->second->u_substep.at(K2)[2][i];
            }
        }

        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u_substep.at(K3)[0][i] = it->second->u[D_UA][i];
                it->second->u_substep.at(K3)[1][i] = it->second->u[D_VA][i];
                it->second->u_substep.at(K3)[2][i] = it->second->u[D_H][i];

                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + dt*it->second->u_substep.at(K3)[0][i];
                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + dt*it->second->u_substep.at(K3)[1][i];
                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + dt*it->second->u_substep.at(K3)[2][i];
            }
        }

        this->Timestep();

        for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
            for (int i = 0; i < it->second->number_bf; i++) {
                it->second->u[UA][i] = it->second->u_substep.at(U0)[0][i] + (dt / 6.0)*
                    (it->second->u_substep.at(K1)[0][i] + 2 * it->second->u_substep.at(K2)[0][i] +
                        2 * it->second->u_substep.at(K3)[0][i] + it->second->u[D_UA][i]);

                it->second->u[VA][i] = it->second->u_substep.at(U0)[1][i] + (dt / 6.0)*
                    (it->second->u_substep.at(K1)[1][i] + 2 * it->second->u_substep.at(K2)[1][i] +
                        2 * it->second->u_substep.at(K3)[1][i] + it->second->u[D_VA][i]);

                it->second->u[H][i] = it->second->u_substep.at(U0)[2][i] + (dt / 6.0)*
                    (it->second->u_substep.at(K1)[2][i] + 2 * it->second->u_substep.at(K2)[2][i] +
                        2 * it->second->u_substep.at(K3)[2][i] + it->second->u[D_H][i]);
            }
        }
    }

    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 3; j++) {
                delete[] it->second->u_substep.at(i)[j];
            }
            delete[] it->second->u_substep.at(i);
        }
        it->second->u_substep.clear();
    }
}

void PROBLEM::WriteDataVTK() {
    std::vector<double> cell_data;

    std::string file_name = "data.vtk";
    std::ofstream file(file_name);

    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        it->second->WriteDataVTK(cell_data, ZB);
    }

    file << "CELL_DATA " << cell_data.size() << '\n';
    file << "SCALARS zb float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (auto it = cell_data.begin(); it != cell_data.end(); it++) {
        file << *it << '\n';
    }

    cell_data.clear();

    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        it->second->WriteDataVTK(cell_data, UA);
    }

    file << "SCALARS ua float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (auto it = cell_data.begin(); it != cell_data.end(); it++) {
        file << *it << '\n';
    }

    cell_data.clear();

    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        it->second->WriteDataVTK(cell_data, VA);
    }

    file << "SCALARS va float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (auto it = cell_data.begin(); it != cell_data.end(); it++) {
        file << *it << '\n';
    }

    cell_data.clear();

    for (auto it = this->mesh->elements.begin(); it != this->mesh->elements.end(); it++) {
        it->second->WriteDataVTK(cell_data, H);
    }

    file << "SCALARS h float 1\n";
    file << "LOOKUP_TABLE default\n";

    for (auto it = cell_data.begin(); it != cell_data.end(); it++) {
        file << *it << '\n';
    }

    cell_data.clear();

    file.close();

    std::string file_name_geom = "geometry.vtk";
    std::string file_name_data = "data.vtk";

    std::ifstream file_geom(file_name_geom, std::ios_base::binary);
    std::ifstream file_data(file_name_data, std::ios_base::binary);

    std::string file_name_merge = "mesh_data.vtk";
    std::ofstream file_merge(file_name_merge, std::ios_base::binary);

    file_merge << file_geom.rdbuf() << file_data.rdbuf();
    file_merge.close();
}