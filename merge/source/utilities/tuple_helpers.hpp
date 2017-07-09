#ifndef TUPLE_HELPERS_HPP
#define TUPLE_HELPERS_HPP

#include "../general_definitions.hpp"

namespace Utilities {
	// has_type implementation from
	//https://stackoverflow.com/questions/25958259/how-do-i-find-out-if-a-tuple-contains-a-type
	template <typename T, typename Tuple>
	struct has_type;

	template <typename T>
	struct has_type<T, std::tuple<>> : std::false_type {};

	template <typename T, typename U, typename... Ts>
	struct has_type<T, std::tuple<U, Ts...>> : has_type<T, std::tuple<Ts...>> {};

	template <typename T, typename... Ts>
	struct has_type<T, std::tuple<T, Ts...>> : std::true_type {};

	// get_index returns the first location of a type
	// https://stackoverflow.com/questions/18063451/get-index-of-a-tuple-elements-type
	template <class T, class Tuple>
	struct index;

	template <class T, class... Ts>
	struct index<T, std::tuple<T, Ts...>> {
		static const std::size_t value = 0;
	};

	template <class T, class U, class... Ts>
	struct index<T, std::tuple<U, Ts...>> {
		static const std::size_t value = 1 + index<T, std::tuple<Ts...>>::value;
	};

	// for_each_in_tuple implementation
	//https://stackoverflow.com/questions/16387354/template-tuple-calling-a-function-on-each-element
	template<int... Is>
	struct seq { };

	template<int N, int... Is>
	struct gen_seq : gen_seq<N - 1, N - 1, Is...> { };

	template<int... Is>
	struct gen_seq<0, Is...> : seq<Is...> { };

	template<typename T, typename F, int... Is>
	void for_each(T&& t, F f, seq<Is...>) {
		//ensure that variables get executed in the correct order
		auto __attribute((unused)) unused = { (f(std::get<Is>(t)), 0)... };
	}

	template<typename... Ts, typename F>
	void for_each_in_tuple(std::tuple<Ts...>& t, F f) {
		for_each(t, f, gen_seq<sizeof...(Ts)>());
	}

	// take two tuples tuple<Pack1...> tuple<Pack2...> and join them to
	// make a new tuple tuple<Pack1..., Pack2...>
	template<typename Tup1, typename Tup2>
	struct tuple_join;

	template<typename... Pack1, typename... Pack2>
	struct tuple_join<std::tuple<Pack1...>, std::tuple<Pack2...> >
	{
		typedef std::tuple<Pack1..., Pack2...> type;
	};
}

#endif