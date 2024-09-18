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

#include <mender-update/update_module/v3/update_module.hpp>

#include <cerrno>
#include <fstream>

#include <unistd.h>

#include <filesystem>

#include <client_shared/conf.hpp>
#include <common/events_io.hpp>
#include <common/io.hpp>
#include <common/log.hpp>
#include <common/path.hpp>
#include <mender-update/context.hpp>

namespace mender {
namespace update {
namespace update_module {
namespace v3 {

using namespace std;

namespace error = mender::common::error;
namespace events = mender::common::events;
namespace expected = mender::common::expected;
namespace conf = mender::client_shared::conf;
namespace context = mender::update::context;
namespace io = mender::common::io;
namespace log = mender::common::log;
namespace path = mender::common::path;

namespace fs = std::filesystem;

error::Error CreateDataFile(
	const fs::path &file_tree_path, const string &file_name, const string &data) {
	string fpath = (file_tree_path / file_name).string();
	auto ex_os = io::OpenOfstream(fpath);
	if (!ex_os) {
		return ex_os.error();
	}
	ofstream &os = ex_os.value();
	if (data != "") {
		auto err = io::WriteStringIntoOfstream(os, data);
		return err;
	}
	return error::NoError;
}

error::Error UpdateModule::PrepareFileTreeDeviceParts(const string &path) {
	// make sure all the required data can be gathered first before creating
	// directories and files
	auto ex_provides = ctx_.LoadProvides();
	if (!ex_provides) {
		return ex_provides.error();
	}

	auto ex_device_type = ctx_.GetDeviceType();
	if (!ex_device_type) {
		return ex_device_type.error();
	}

	const fs::path file_tree_path {path};

	const fs::path tmp_subdir_path = file_tree_path / "tmp";
	auto err = path::CreateDirectories(tmp_subdir_path.string());
	if (err != error::NoError) {
		return err;
	}

	auto provides = ex_provides.value();
	auto write_provides_into_file = [&](const string &key) {
		return CreateDataFile(
			file_tree_path,
			"current_" + key,
			(provides.count(key) != 0) ? provides[key] + "\n" : "");
	};

	err = CreateDataFile(file_tree_path, "version", "3\n");
	if (err != error::NoError) {
		return err;
	}

	err = write_provides_into_file("artifact_name");
	if (err != error::NoError) {
		return err;
	}
	err = write_provides_into_file("artifact_group");
	if (err != error::NoError) {
		return err;
	}

	auto device_type = ex_device_type.value();
	err = CreateDataFile(file_tree_path, "current_device_type", device_type + "\n");
	if (err != error::NoError) {
		return err;
	}

	return error::NoError;
}

error::Error UpdateModule::CleanAndPrepareFileTree(
	const string &path, artifact::PayloadHeaderView &payload_meta_data) {
	const fs::path file_tree_path {path};

	std::error_code ec;
	fs::remove_all(file_tree_path, ec);
	if (ec) {
		return error::Error(
			ec.default_error_condition(), "Could not clean File Tree for Update Module");
	}

	auto err = PrepareFileTreeDeviceParts(path);
	if (err != error::NoError) {
		return err;
	}

	//
	// Header
	//

	const fs::path header_subdir_path = file_tree_path / "header";
	err = path::CreateDirectories(header_subdir_path.string());
	if (err != error::NoError) {
		return err;
	}

	err = CreateDataFile(
		header_subdir_path, "artifact_group", payload_meta_data.header.artifact_group);
	if (err != error::NoError) {
		return err;
	}

	err =
		CreateDataFile(header_subdir_path, "artifact_name", payload_meta_data.header.artifact_name);
	if (err != error::NoError) {
		return err;
	}

	err = CreateDataFile(header_subdir_path, "payload_type", payload_meta_data.header.payload_type);
	if (err != error::NoError) {
		return err;
	}

	err = CreateDataFile(
		header_subdir_path, "header-info", payload_meta_data.header.header_info.verbatim.Dump());
	if (err != error::NoError) {
		return err;
	}

	err = CreateDataFile(
		header_subdir_path, "type-info", payload_meta_data.header.type_info.verbatim.Dump());
	if (err != error::NoError) {
		return err;
	}

	err =
		CreateDataFile(header_subdir_path, "meta-data", payload_meta_data.header.meta_data.Dump());
	if (err != error::NoError) {
		return err;
	}

	// Make sure all changes are permanent, even across spontaneous reboots. We don't want to
	// have half a tree when trying to recover from that.
	return path::DataSyncRecursively(path);
}

error::Error UpdateModule::EnsureRootfsImageFileTree(const string &path) {
	// Historical note: Versions of the client prior to 4.0 had the rootfs-image module built
	// in. Because of this it has no Update Module File Tree. So if we are upgrading, we might
	// hit an on-going upgrade without a File Tree. It's too late to create a complete one with
	// all the artifact content by the time we get here, but at least we can create one which
	// has the current Provides information, as well as a folder for the Update Module to run
	// in.
	ifstream payload_type(path::Join(path, "header", "payload_type"));
	if (payload_type.good()) {
		string type;
		payload_type >> type;
		if (payload_type.good() && type == "rootfs-image") {
			// If we have a File Tree with the rootfs-image type, we assume we are
			// fine. This is actually not completely safe in an upgrade situation,
			// because the old <4.0 client will not have cleaned the tree, and it could
			// be old. However, this will *only* happen in an upgrade situation from
			// <4.0 to >=4.0, and I can't think of a way it could be exploited. Also,
			// the rootfs-image module does not use any of the information ATM.
			return error::NoError;
		}
	}

	return PrepareFileTreeDeviceParts(path);
}

error::Error UpdateModule::DeleteFileTree(const string &path) {
	try {
		fs::remove_all(fs::path {path});
	} catch (const fs::filesystem_error &e) {
		return error::Error(
			e.code().default_error_condition(),
			"Failed to recursively remove directory '" + path + "': " + e.what());
	}

	return error::NoError;
}

expected::ExpectedStringVector DiscoverUpdateModules(const conf::MenderConfig &config) {
	vector<string> ret {};
	fs::path file_tree_path = fs::path(config.paths.GetDataStore()) / "modules/v3";

	try {
		for (auto entry : fs::directory_iterator(file_tree_path)) {
			const fs::path file_path = entry.path();
			const string file_path_str = file_path.string();
			if (!fs::is_regular_file(file_path)) {
				log::Warning("'" + file_path_str + "' is not a regular file");
				continue;
			}

			const fs::perms perms = entry.status().permissions();
			if ((perms & (fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec))
				== fs::perms::none) {
				log::Warning("'" + file_path_str + "' is not executable");
				continue;
			}

			ret.push_back(file_path_str);
		}
	} catch (const fs::filesystem_error &e) {
		auto code = e.code();
		if (code.value() == ENOENT) {
			// directory not found is not an error, just return an empty vector
			return ret;
		}
		// everything (?) else is an error
		return expected::unexpected(error::Error(
			code.default_error_condition(),
			"Failed to discover update modules in '" + file_tree_path.string() + "': " + e.what()));
	}

	return ret;
}

error::Error UpdateModule::PrepareStreamNextPipe() {
	download_->stream_next_path_ = path::Join(update_module_workdir_, "stream-next");

	if (::mkfifo(download_->stream_next_path_.c_str(), 0600) != 0) {
		int err = errno;
		return error::Error(
			generic_category().default_error_condition(err),
			"Unable to create `stream-next` at " + download_->stream_next_path_);
	}
	return error::NoError;
}

error::Error UpdateModule::OpenStreamNextPipe(ExpectedWriterHandler open_handler) {
	auto opener = make_shared<AsyncFifoOpener>(download_->event_loop_);
	download_->stream_next_opener_ = opener;
	return opener->AsyncOpen(download_->stream_next_path_, open_handler);
}

error::Error UpdateModule::PrepareAndOpenStreamPipe(
	const string &path, ExpectedWriterHandler open_handler) {
	auto fs_path = fs::path(path);
	std::error_code ec;
	if (!fs::create_directories(fs_path.parent_path(), ec) && ec) {
		return error::Error(
			ec.default_error_condition(),
			"Could not create stream directory at " + fs_path.parent_path().string());
	}

	if (::mkfifo(path.c_str(), 0600) != 0) {
		int err = errno;
		return error::Error(
			generic_category().default_error_condition(err),
			"Could not create stream FIFO at " + path);
	}

	auto opener = make_shared<AsyncFifoOpener>(download_->event_loop_);
	download_->current_stream_opener_ = opener;
	return opener->AsyncOpen(path, open_handler);
}

error::Error UpdateModule::PrepareDownloadDirectory(const string &path) {
	auto fs_path = fs::path(path);
	std::error_code ec;
	if (!fs::create_directories(fs_path, ec) && ec) {
		return error::Error(
			ec.default_error_condition(), "Could not create `files` directory at " + path);
	}

	return error::NoError;
}

error::Error UpdateModule::DeleteStreamsFiles() {
	try {
		fs::path p {download_->stream_next_path_};
		fs::remove_all(p);

		p = fs::path(update_module_workdir_) / "streams";
		fs::remove_all(p);
	} catch (fs::filesystem_error &e) {
		return error::Error(
			e.code().default_error_condition(), "Could not remove " + download_->stream_next_path_);
	}

	return error::NoError;
}

AsyncFifoOpener::AsyncFifoOpener(events::EventLoop &loop) :
	event_loop_ {loop},
	cancelled_ {make_shared<bool>(true)},
	destroying_ {make_shared<bool>(false)} {
}

AsyncFifoOpener::~AsyncFifoOpener() {
	*destroying_ = true;
	Cancel();
}

error::Error AsyncFifoOpener::AsyncOpen(const string &path, ExpectedWriterHandler handler) {
	// Excerpt from fifo(7) man page:
	// ------------------------------
	// The FIFO must be opened on both ends (reading and writing) before data can be passed.
	// Normally, opening the FIFO blocks until the other end is opened also.
	//
	// A process can open a FIFO in nonblocking mode. In this case, opening for read-only
	// succeeds even if no one has opened on the write side yet and opening for write-only fails
	// with ENXIO (no such device or address) unless the other end has already been opened.
	//
	// Under Linux, opening a FIFO for read and write will succeed both in blocking and
	// nonblocking mode. POSIX leaves this behavior undefined.  This can be used to open a FIFO
	// for writing while there are no readers available.
	// ------------------------------
	//
	// We want to open the pipe to check if a process is going to read from it. But we cannot do
	// this in a blocking fashion, because we are also waiting for the process to terminate,
	// which happens for Update Modules that want the client to download the files for them. So
	// we need this AsyncFifoOpener class, which does the work in a thread.

	if (!*cancelled_) {
		return error::Error(
			make_error_condition(errc::operation_in_progress), "Already running AsyncFifoOpener");
	}

	*cancelled_ = false;
	path_ = path;
	thread_ = thread([this, handler]() {
		auto writer = make_shared<events::io::AsyncFileDescriptorWriter>(event_loop_);
		// This will block for as long as there are no FIFO readers.
		auto err = writer->Open(path_);

		io::ExpectedAsyncWriterPtr exp_writer;
		if (err != error::NoError) {
			exp_writer = expected::unexpected(err);
		} else {
			exp_writer = writer;
		}

		auto &cancelled = cancelled_;
		auto &destroying = destroying_;
		event_loop_.Post([handler, exp_writer, cancelled, destroying]() {
			if (*destroying) {
				return;
			}

			if (*cancelled) {
				handler(expected::unexpected(error::Error(
					make_error_condition(errc::operation_canceled), "AsyncFifoOpener cancelled")));
				return;
			}

			handler(exp_writer);
		});
	});

	return error::NoError;
}

void AsyncFifoOpener::Cancel() {
	if (*cancelled_) {
		return;
	}

	*cancelled_ = true;

	// Open the other end of the pipe to jerk the first end loose.
	int fd = ::open(path_.c_str(), O_RDONLY | O_NONBLOCK);
	if (fd < 0) {
		int errnum = errno;
		log::Error(string("Cancel::open() returned error: ") + strerror(errnum));
	}

	thread_.join();

	::close(fd);
}

} // namespace v3
} // namespace update_module
} // namespace update
} // namespace mender
