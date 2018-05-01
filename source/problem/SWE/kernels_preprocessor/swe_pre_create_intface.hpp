#ifndef SWE_PRE_CREATE_INTFACE_HPP
#define SWE_PRE_CREATE_INTFACE_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_interfaces_kernel(
    ProblemMeshType& mesh,
    std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
    Writer<SWE::Problem>& writer) {
    // *** //

    using InterfaceTypes = Geometry::InterfaceTypeTuple<SWE::Data, SWE::IS::Interface, SWE::IS::Levee>;

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); it++) {
        switch (it->first) {
            case SWE::BoundaryConditions::internal:
                using InterfaceTypeInterface = std::tuple_element<0, InterfaceTypes>::type;

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                    if (it->second.find(key_pre_int_ex) != it->second.end()) {
                        mesh.template CreateInterface<InterfaceTypeInterface>(itt->second,
                                                                              it->second.find(key_pre_int_ex)->second);
                    }

                    it->second.erase(itt);
                }

                break;
            case SWE::BoundaryConditions::internal_barrier:
                using InterfaceTypeLevee = std::tuple_element<1, InterfaceTypes>::type;

                for (auto itt = it->second.begin(); itt != it->second.end(); itt++) {
                    std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                    if (it->second.find(key_pre_int_ex) != it->second.end()) {
                        mesh.template CreateInterface<InterfaceTypeLevee>(itt->second,
                                                                          it->second.find(key_pre_int_ex)->second);
                    }

                    it->second.erase(itt);
                }

                break;
        }

        if (writer.WritingLog()) {
            writer.GetLogFile() << "Number of interfaces: " << mesh.GetNumberInterfaces() << std::endl;
        }
    }
}
}

#endif