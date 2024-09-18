// Copyright 2023 Northern.tech AS
//
//    Licensed under the Apache License, Version 2.0 (the "License");
//    you may not use this file except in compliance with the License.
//    You may obtain a copy of the License at
//
//        http://www.apache.org/licenses/LICENSE-2.0
//
//    Unless required by applicable law or agreed to in writing, software
//    distributed under the License is distributed on an "AS IS" BASIS,
//    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
//    See the License for the specific language governing permissions and
//    limitations under the License.

#ifndef MENDER_COMMON_EXPECTED_HPP
#define MENDER_COMMON_EXPECTED_HPP

#ifdef __cpp_lib_expected
#include <expected>
#else
#include <tl/expected.hpp>
#endif

#include <common/error.hpp>

#include <cassert>
#include <cstdint>
#include <string>
#include <vector>
#include <unordered_set>


namespace mender {
namespace common {
namespace expected {

using namespace std;

#ifdef __cpp_lib_expected
using std::expected;
using std::unexpected;
#else
using tl::expected;
template <typename V>
tl::unexpected<V> unexpected(V &&v) {
	return tl::make_unexpected(v);
}
template <typename V>
tl::unexpected<V> unexpected(V &v) {
	return tl::make_unexpected(v);
}
#endif

template <typename T>
using Expected = expected<T, error::Error>;

using ExpectedString = Expected<string>;
using ExpectedBytes = Expected<vector<uint8_t>>;
using ExpectedInt = Expected<int>;
using ExpectedInt64 = Expected<int64_t>;
using ExpectedDouble = Expected<double>;
using ExpectedBool = Expected<bool>;
using ExpectedSize = Expected<size_t>;
using ExpectedBool = Expected<bool>;
using ExpectedLong = Expected<long>;
using ExpectedLongLong = Expected<long long>;
using ExpectedStringVector = Expected<vector<string>>;

template <typename T>
using ExpectedVector = Expected<vector<T>>;

template <typename T>
using ExpectedUnorderedSet = Expected<unordered_set<T>>;

} // namespace expected
} // namespace common
} // namespace mender

#endif // MENDER_COMMON_EXPECTED_HPP
