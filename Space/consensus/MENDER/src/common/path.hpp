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

#ifndef MENDER_COMMON_PATH_HPP
#define MENDER_COMMON_PATH_HPP

#include <functional>
#include <string>

#include <common/error.hpp>
#include <common/expected.hpp>
#include <common/log.hpp>

namespace mender {
namespace common {
namespace path {

using namespace std;

namespace error = mender::common::error;
namespace expected = mender::common::expected;

enum class Perms {
	Owner_read,
	Owner_write,
	Owner_exec,
	Group_read,
	Group_write,
	Group_exec,
	Others_read,
	Others_write,
	Others_exec,
};

error::Error Permissions(const string &file_path, vector<Perms> perms);

string JoinOne(const string &prefix, const string &path);

template <typename... Paths>
string Join(const string &prefix, const Paths &...paths) {
	string final_path {prefix};
	for (const auto &path : {paths...}) {
		final_path = JoinOne(final_path, path);
	}
	return final_path;
}

string BaseName(const string &path);
string DirName(const string &path);

expected::ExpectedString Canonical(const string &path);

bool IsAbsolute(const string &path);

bool FileExists(const string &path);

error::Error FileDelete(const string &path);

error::Error DeleteRecursively(const string &path);

expected::ExpectedBool IsExecutable(const string &path, const bool warn = true);

expected::ExpectedBool AreFilesIdentical(const string &file_one, const string &file_two);

expected::ExpectedUnorderedSet<string> ListFiles(const string &in, function<bool(string)> matcher);

error::Error CreateDirectory(const string &path);

error::Error CreateDirectories(const string &dir);

error::Error DataSyncRecursively(const string &dir);

error::Error Rename(const string &oldname, const string &newname);
error::Error FileCopy(const string &what, const string &where);

} // namespace path
} // namespace common
} // namespace mender

#endif // MENDER_COMMON_PATH_HPP
