#ifndef EHDG_SWE_PRE_CREATE_INTFACE_HPP
#define EHDG_SWE_PRE_CREATE_INTFACE_HPP

namespace SWE {
namespace EHDG {
template <typename RawBoundaryType>
void Problem::create_interfaces(std::map<uchar, std::map<std::pair<uint, uint>, RawBoundaryType>>& raw_boundaries,
                                ProblemMeshType& mesh,
                                ProblemInputType& problem_input,
                                Writer<Problem>& writer) {
    // *** //
    using InterfaceTypes = Geometry::InterfaceTypeTuple<Data, IS::Internal>;

    for (auto it = raw_boundaries.begin(); it != raw_boundaries.end(); ++it) {
        if (it->first == SWE::BoundaryTypes::internal) {
            using InterfaceTypeInternal = std::tuple_element<0, InterfaceTypes>::type;

            uint n_intface_old_internal = mesh.GetNumberInterfaces();

            auto itt = it->second.begin();
            while (itt != it->second.end()) {
                std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{itt->first.second, itt->first.first};

                if (it->second.find(key_pre_int_ex) != it->second.end()) {
                    auto& raw_boundary_in = itt->second;
                    auto& raw_boundary_ex = it->second.find(key_pre_int_ex)->second;

                    mesh.template CreateInterface<InterfaceTypeInternal>(std::move(raw_boundary_in),
                                                                         std::move(raw_boundary_ex));
                }

                it->second.erase(itt++);
            }

            if (writer.WritingLog()) {
                writer.GetLogFile() << "Number of internal interfaces: "
                                    << mesh.GetNumberInterfaces() - n_intface_old_internal << std::endl;
            }
        }
    }

    mesh.CallForEachInterface([](auto& intface) { intface.specialization.Initialize(intface); });
}
}
}

#endif
