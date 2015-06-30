#ifndef GWIZ_VCFFILEREADER_H
#define GWIZ_VCFFILEREADER_H

#include "core/variant/IVariantList.h"
#include "core/file/IFile.h"

#include "core/parser/VCFParser.hpp"
#include "core/parser/ChromParser.hpp"
#include "core/region/Region.h"
#include "core/util/SharedCreator.hpp"

#include "Variant.h"

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>

#include <list>
#include <tuple>
#include <map>
#include <thread>
#include <mutex>
#include <future>

namespace gwiz
{
	class VCFFileReader : private boost::noncopyable, public SharedCreator< VCFFileReader >
	{
    public:
		typedef std::shared_ptr<VCFFileReader> SharedPtr;
		typedef std::weak_ptr<VCFFileReader> WeakPtr;

		inline static VCFFileReader::SharedPtr CreateVCFFileReader(const std::string& path)
		{
			/*
			_private_const p;
			_private_destr d;
			auto vcfFileReaderPtr = std::make_shared< VCFFileReader >(p, path, d);
			vcfFileReaderPtr->m_this_wk_ptr = vcfFileReaderPtr;
			return vcfFileReaderPtr;
			*/
			return nullptr;
		}

		uint32_t getID() { return this->m_id; }
		std::vector< IVariant::SharedPtr > getVariantsInRegion(Region::SharedPtr regionPtr);
	protected:
		VCFFileReader(const std::string& path);
		~VCFFileReader();
	private:

		void Open();
        void readHeader();
		position getPositionFromLine(const char* line);
		void setFileReader(const std::string& path);

		VCFFileReader::WeakPtr m_this_wk_ptr;
		IFile::SharedPtr m_file_ptr;
		VariantParser< const char* > m_vcf_parser;
		uint32_t m_id;

		std::mutex m_region_mutex;

	};
} // end namespace gwiz

#endif  //GWIZ_VCFFILEREADER_H
