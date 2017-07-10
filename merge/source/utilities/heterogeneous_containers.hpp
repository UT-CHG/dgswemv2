#ifndef HETEROGENEOUS_CONTAINERS_HPP
#define HETEROGENEOUS_CONTAINERS_HPP

#include "tuple_helpers.hpp"

namespace Utilities {
	template<typename... Ts>
	struct HeterogeneousVector {
		using TupleType = std::tuple<Ts...>;
		std::tuple<std::vector<Ts>...> data;

		uint size() {
			uint size = 0;

			for_each_in_tuple(this->data, [&size](const auto& vector) { size += vector.size(); });

			return size;
		}

		template<typename T>
		void emplace_back(T&& t) {
			static_assert(has_type<T, TupleType>::value,
				"Error in HeterogeneousVector::emplace_back: Type not found");

			std::get<index<T, TupleType>::value>(data).emplace_back(std::forward<T>(t));
		}

		template<typename T>
		T& at(uint i) {
			static_assert(has_type<T, TupleType>::value,
				"Error in HeterogeneousVector::at: Type not found");

			return std::get<index<T, TupleType>::value>(data).at(i);
		}
	};

	template<typename... Ts>
	struct HeterogeneousMap {
		using TupleType = std::tuple<Ts...>;
		std::tuple<std::map<uint, Ts>...> data;

		uint size() {
			uint size = 0;

			for_each_in_tuple(this->data, [&size](const auto& map) { size += map.size(); });

			return size;
		}

		template<typename T>
		void emplace(uint n, T&& t) {
			static_assert(has_type<T, TupleType>::value,
				"Error in HeterogeneousMap::emplace: Type not found");

			std::get<index<T, TupleType>::value>(data).emplace(std::make_pair(n, std::forward<T>(t)));
		}

		template<typename T>
		T& at(uint i) {
			static_assert(has_type<T, TupleType>::value,
				"Error in HeterogeneousMap::at: Type not found");

			return std::get<index<T, TupleType>::value>(data).at(i);
		}
	};
}

#endif