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

#include <common/path.hpp>

#include <filesystem>
#include <string>
#include <unordered_set>

#include <common/error.hpp>

namespace mender {
namespace common {
namespace path {

using namespace std;
namespace fs = std::filesystem;

unordered_map<Perms, fs::perms> perm_map = {
	{Perms::Owner_exec, fs::perms::owner_exec},
	{Perms::Owner_read, fs::perms::owner_read},
	{Perms::Owner_write, fs::perms::owner_write},
	{Perms::Group_read, fs::perms::group_read},
	{Perms::Group_write, fs::perms::group_write},
	{Perms::Group_exec, fs::perms::group_exec},
	{Perms::Others_read, fs::perms::others_read},
	{Perms::Others_write, fs::perms::others_write},
	{Perms::Others_exec, fs::perms::others_exec},
};

string JoinOne(const string &prefix, const string &suffix) {
	return (fs::path(prefix) / suffix).string();
}

string BaseName(const string &path) {
	return fs::path(path).filename().string();
}

string DirName(const string &path) {
	return fs::path(path).parent_path().string();
}

expected::ExpectedString Canonical(const string &path) {
	error_code ec;
	auto canonical = fs::canonical(path, ec);
	if (ec) {
		return expected::unexpected(
			error::Error(ec.default_error_condition(), "Could not get canonical path"));
	}
	return canonical.string();
}

bool IsAbsolute(const string &path) {
	return fs::path(path).is_absolute();
}

bool FileExists(const string &path) {
	try {
		return fs::exists(path);
	} catch (const fs::filesystem_error &e) {
		log::Error("Could not check file existence of '" + path + "': " + e.what());
		return false;
	}
}

error::Error FileDelete(const string &path) {
	error_code ec;
	bool deleted = fs::remove(fs::path(path), ec);
	if (not deleted) {
		return error::Error(
			ec.default_error_condition(),
			"Failed to remove the file: '" + path + "'. error: " + ec.message());
	}
	return error::NoError;
}

error::Error DeleteRecursively(const string &path) {
	error_code ec;
	fs::remove_all(path, ec);
	if (ec) {
		return error::Error(ec.default_error_condition(), "Could not remove path");
	}
	return error::NoError;
}

expected::ExpectedBool IsExecutable(const string &file_path, const bool warn) {
	try {
		fs::perms perms = fs::status(file_path).permissions();
		if ((perms & (fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec))
			== fs::perms::none) {
			if (warn) {
				log::Warning("'" + file_path + "' is not executable");
			}
			return false;
		}
		return true;
	} catch (const fs::filesystem_error &e) {
		return expected::unexpected(error::Error(
			e.code().default_error_condition(),
			"Could not check executable status of '" + file_path + "'"));
	}
}

error::Error Permissions(const string &file_path, const vector<Perms> perms) {
	if (perms.size() == 0) {
		return error::NoError;
	}
	fs::perms p {};
	std::for_each(perms.cbegin(), perms.cend(), [&p](const Perms perm) { p |= perm_map.at(perm); });
	try {
		fs::permissions(file_path, p);
	} catch (const fs::filesystem_error &e) {
		return error::Error(
			e.code().default_error_condition(), "Could not set permissions on '" + file_path + "'");
	}
	return error::NoError;
}

expected::ExpectedUnorderedSet<string> ListFiles(
	const string &in_directory, function<bool(string)> matcher) {
	try {
		unordered_set<string> matching_files {};
		fs::path dir_path(in_directory);
		if (!fs::exists(dir_path)) {
			auto err {errno};
			return expected::unexpected(error::Error(
				generic_category().default_error_condition(err),
				"No such file or directory: " + in_directory));
		}

		for (const auto &entry : fs::directory_iterator {dir_path}) {
			fs::path file_path = entry.path();
			if (!fs::is_regular_file(file_path)) {
				log::Warning("'" + file_path.string() + "'" + " is not a regular file. Ignoring.");
				continue;
			}

			if (matcher(file_path)) {
				matching_files.insert(file_path);
			}
		}

		return matching_files;
	} catch (const fs::filesystem_error &e) {
		return expected::unexpected(error::Error(
			e.code().default_error_condition(), "Could not list files in '" + in_directory + "'"));
	}
}

error::Error CreateDirectory(const string &path) {
	try {
		fs::path fs_path {path};
		if (not fs::create_directory(fs_path)) {
			auto err {errno};
			return error::Error(
				generic_category().default_error_condition(err),
				"Failed to create the directory: " + path);
		}
	} catch (const fs::filesystem_error &e) {
		return error::Error(
			e.code().default_error_condition(), "Failed to create directory: '" + path + "'");
	}
	return error::NoError;
}

error::Error CreateDirectories(const string &dir) {
	try {
		const fs::path p {dir};
		fs::create_directories(p);
	} catch (const fs::filesystem_error &e) {
		return error::Error(
			e.code().default_error_condition(), "Failed to create directory: '" + dir + "'");
	}
	return error::NoError;
}

error::Error Rename(const string &oldname, const string &newname) {
	error_code ec;
	fs::rename(fs::path(oldname), fs::path(newname), ec);
	if (ec) {
		return error::Error(
			ec.default_error_condition(),
			"Failed to rename '" + oldname + "' to '" + newname + "'. error: " + ec.message());
	}
	return error::NoError;
}

error::Error FileCopy(const string &what, const string &where) {
	error_code ec;
	fs::copy_file(fs::path(what), fs::path(where), fs::copy_options::overwrite_existing, ec);
	if (ec) {
		return error::Error(
			ec.default_error_condition(),
			"Failed to copy '" + what + "' to '" + where + "'. error: " + ec.message());
	}
	return error::NoError;
}

} // namespace path
} // namespace common
} // namespace mender
