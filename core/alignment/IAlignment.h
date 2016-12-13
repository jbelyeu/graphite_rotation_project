#ifndef GRAPHITE_IALIGNMENT_H
#define GRAPHITE_IALIGNMENT_H

#include "core/util/Types.h"
#include "core/util/Noncopyable.hpp"
#include "core/allele/IAllele.h"

#include <memory>
#include <unordered_map>
#include <mutex>
#include <vector>
#include <string>

namespace graphite
{

	class IMapping;
	class Sample;
	class IAlignment : private Noncopyable
	{
	public:
		typedef std::shared_ptr< IAlignment > SharedPtr;

	    IAlignment() : m_mapping_mutex(new std::recursive_mutex()) {}
		virtual ~IAlignment() { delete this->m_mapping_mutex; }

		virtual const char* getSequence() = 0;
		virtual const position getPosition() = 0;
		virtual const size_t getLength() = 0;
		virtual const uint32_t getFlag() {return this->m_flag; };
		virtual const uint32_t getTemplateLength() {return this->m_template_len; };
		virtual const uint32_t getQual() {return this->m_qual; };
		virtual const uint32_t getChrID() {return this->m_chr_id; };
		virtual const uint32_t getPositionNext() {return this->m_pos_next; };
		virtual const uint32_t getRefNext() {return this->m_ref_next; };
		virtual const std::string getID() { return ""; }
		virtual const bool isFirstMate() { return false;}
		virtual const bool isMapped() { return false; }
		virtual const bool isReverseStrand() { return false; }
		virtual const bool isDuplicate() { return false; }
		virtual const uint16_t getOriginalMapQuality() { return 0; }
		virtual std::weak_ptr< IMapping > getMapping() { return this->m_mapping_wptr; }
		virtual void setMapping(std::shared_ptr< IMapping > mappingPtr)
		{
			std::lock_guard< std::recursive_mutex > r_lock(*this->m_mapping_mutex);
			this->m_mapping_wptr = mappingPtr;
		}
		std::recursive_mutex* getMappingMutex() { return this->m_mapping_mutex; }
		const std::shared_ptr< Sample > getSample() { return m_sample_ptr; }

		virtual const void setSequence(char* seq, uint32_t len) = 0;
		virtual const void setFlag(uint32_t flag) { m_flag = flag; }
		virtual const void setTemplateLength(uint32_t tlen) { m_template_len = tlen; }
		virtual const void setQual(uint32_t qual) { m_qual = qual; }
		virtual const void setChrID(uint32_t chr_id) { m_chr_id = chr_id; }
		virtual const void setPositionNext(uint32_t pos_next) { m_pos_next = pos_next; }
		virtual const void setRefNext(uint32_t ref_next) { m_ref_next = ref_next; }
		virtual const void removeSequence() = 0;
		virtual const void incrementReferenceCount() = 0;

	protected:
		std::mutex m_mutex;
		std::unordered_map< uint32_t, std::string > m_mapped_variants_information;
		std::weak_ptr< IMapping > m_mapping_wptr;
		std::recursive_mutex* m_mapping_mutex;
		std::shared_ptr< Sample > m_sample_ptr;
		uint32_t m_flag;
		uint32_t m_template_len;
		uint32_t m_qual;
		uint32_t m_chr_id;
		uint32_t m_pos_next;
		uint32_t m_ref_next;
	};
}

#endif //GRAPHITE_IALIGNMENT_H
