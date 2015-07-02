#ifndef GWIZ_ALLELE_H
#define GWIZ_ALLELE_H

#include "core/sequence/Sequence.h"
#include "IAllele.h"
#include "AlleleMetaData.h"

namespace gwiz
{
	class VCFFileReader;
	class Allele : public IAllele
	{
	public:
		typedef std::shared_ptr< Allele > SharedPtr;
	    Allele(std::shared_ptr< Sequence > sequencePtr) :
		    m_sequence_ptr(sequencePtr)
		{
		}

		~Allele() {}

		IAllele::SharedPtr copyAllele() override
		{
			auto allelePtr = std::make_shared< Allele >(this->m_sequence_ptr);
			allelePtr->m_allele_meta_data_ptr = this->m_allele_meta_data_ptr;
			return allelePtr;
		}

		std::shared_ptr< Sequence > getSequencePtr() override { return this->m_sequence_ptr; }
		const char* getSequence() override { return this->m_sequence_ptr->getSequence(); }
		std::string getSequenceString() override { return this->m_sequence_ptr->getSequenceString(); }
		void setMetaData(AlleleMetaData::SharedPtr alleleMetaDataPtr) { this->m_allele_meta_data_ptr = alleleMetaDataPtr; }
		AlleleMetaData::SharedPtr getMetaData() { return this->m_allele_meta_data_ptr; }

	protected:
		Allele() {}
		std::shared_ptr< Sequence > m_sequence_ptr;
		AlleleMetaData::SharedPtr m_allele_meta_data_ptr;
	};
}

#endif //GWIZ_ALLELE_H
