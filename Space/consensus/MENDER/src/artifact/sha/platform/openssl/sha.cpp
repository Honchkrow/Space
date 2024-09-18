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

#include <cstdint>
#include <string>
#include <ctime>
#include <iomanip>
#include <vector>
#include <sstream>
#include <algorithm>

#include <openssl/evp.h>
#include <artifact/sha/sha.hpp>

#include <common/common.hpp>
#include <common/io.hpp>


namespace mender {
namespace sha {

static const size_t SHA_256_digest_length = 32;

namespace log = mender::common::log;
namespace io = mender::common::io;


Reader::Reader(io::Reader &reader, const std::string &expected_sha) :
	sha_handle_(EVP_MD_CTX_new(), [](EVP_MD_CTX *ctx) { EVP_MD_CTX_free(ctx); }),
	wrapped_reader_ {reader},
	expected_sha_ {expected_sha} {
	if (EVP_DigestInit_ex(sha_handle_.get(), EVP_sha256(), nullptr) != 1) {
		log::Error("Failed to initialize the shasummer");
		initialized_ = false;
		return;
	}
	initialized_ = true;
}

expected::ExpectedSize Reader::Read(
	vector<uint8_t>::iterator start, vector<uint8_t>::iterator end) {
	if (!initialized_) {
		return expected::unexpected(MakeError(
			InitializationError,
			"The ShaReader was not properly initialized. Shasumming is not possible"));
	}

	auto bytes_read = wrapped_reader_.Read(start, end);
	if (!bytes_read) {
		return bytes_read;
	}

	// bytes_read == 0 == EOF marker in our Reader/Writer interface implementation
	if (bytes_read.value() == 0) {
		auto real_sha = this->ShaSum();
		if (!real_sha) {
			return expected::unexpected(real_sha.error());
		}
		if (expected_sha_.size() > 0 and real_sha.value() != expected_sha_) {
			return expected::unexpected(MakeError(
				ShasumMismatchError,
				"The checksum of the read byte-stream does not match the expected checksum, (expected): "
					+ expected_sha_ + " (calculated): " + real_sha.value().String()));
		}
		this->done_ = true;
		this->shasum_ = real_sha.value();
		return 0;
	}

	if (EVP_DigestUpdate(sha_handle_.get(), &start[0], bytes_read.value()) != 1) {
		return expected::unexpected(MakeError(ShasumCreationError, "Failed to create the shasum"));
	}

	return bytes_read.value();
}


ExpectedSHA Reader::ShaSum() {
	if (!initialized_) {
		return expected::unexpected(MakeError(
			InitializationError,
			"The ShaReader was not properly initialized. Shasumming is not possible"));
	}
	if (done_) {
		return this->shasum_;
	}

	vector<uint8_t> hash(EVP_MAX_MD_SIZE);
	unsigned int hash_length = 0;

	if (EVP_DigestFinal_ex(sha_handle_.get(), hash.data(), &hash_length) != 1) {
		return expected::unexpected(
			MakeError(ShasumCreationError, "Failed to create the shasum. OpenSSL error: "));
	}

	if (hash_length != SHA_256_digest_length) {
		return expected::unexpected(MakeError(
			ShasumCreationError,
			"SHA of unexpected length: " + std::to_string(hash_length) + " expected length: 32"));
	}

	return SHA(hash, SHA_256_digest_length);
}

} // namespace sha
} // namespace mender
