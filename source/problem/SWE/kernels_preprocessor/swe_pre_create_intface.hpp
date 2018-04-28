#ifndef SWE_PRE_CREATE_INTFACE_HPP
#define SWE_PRE_CREATE_INTFACE_HPP

namespace SWE {
template <typename RawBoundaryType>
void Problem::create_interfaces_kernel(ProblemMeshType&                                  mesh,
                                       std::map<std::pair<uint, uint>, RawBoundaryType>& pre_interfaces,
                                       Writer<SWE::Problem>&                             writer) {
    using InterfaceType = std::tuple_element<0, Geometry::InterfaceTypeTuple<SWE::Data, SWE::IS::Empty>>::type;

    for (auto it = pre_interfaces.begin(); it != pre_interfaces.end(); it++) {
        std::pair<uint, uint> key_pre_int_ex = std::pair<uint, uint>{it->first.second, it->first.first};

        if (pre_interfaces.find(key_pre_int_ex) != pre_interfaces.end()) {
            mesh.template CreateInterface<InterfaceType>(it->second, pre_interfaces.find(key_pre_int_ex)->second);
        }

        pre_interfaces.erase(it);
    }

    if (writer.WritingLog()) {
        writer.GetLogFile() << "Number of interfaces: " << mesh.GetNumberInterfaces() << std::endl;
    }
}
}

#endif