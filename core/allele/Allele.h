#ifndef GRAPHITE_ALLELE_H
#define GRAPHITE_ALLELE_H

#include "core/sequence/SequenceManager.h"
#include "IAllele.h"
#include "AlleleMetaData.h"
#include "core/alignment/IAlignment.h"
#include "core/alignment/Sample.hpp"

namespace graphite
{
	class VCFFileReader;
	class IAlignmentReader;
	class Allele : public IAllele
	{
	public:
		typedef std::shared_ptr< Allele > SharedPtr;
	    Allele(Sequence::SharedPtr sequencePtr) :
		    m_sequence_ptr(sequencePtr),
			m_allele_meta_data_ptr(std::make_shared< AlleleMetaData >(0, 0))
		{
		}

	    Allele(const std::string& seqString) :
			m_sequence_ptr(SequenceManager::Instance()->getSequence(seqString)),
			m_allele_meta_data_ptr(std::make_shared< AlleleMetaData >(0, 0))
		{
		}

	    Allele(const std::string& seqString, AlleleMetaData::SharedPtr alleleMetaDataPtr) :
			m_sequence_ptr(SequenceManager::Instance()->getSequence(seqString)),
			m_allele_meta_data_ptr(alleleMetaDataPtr)
		{
		}

		~Allele()
		{
		}


		size_t getLength() override { return this->m_sequence_ptr->getLength(); }
		std::shared_ptr< Sequence > getSequencePtr() override { return this->m_sequence_ptr; }
		const char* getSequence() override { return this->m_sequence_ptr->getSequence(); }
		void setSequence(std::shared_ptr< Sequence > sequencePtr) override { this->m_sequence_ptr = sequencePtr; }
		std::string getSequenceString() override { return this->m_sequence_ptr->getSequenceString(); }
		virtual void setAlleleMetaData(AlleleMetaData::SharedPtr alleleMetaDataPtr)  override { this->m_allele_meta_data_ptr = alleleMetaDataPtr; }
		virtual AlleleMetaData::SharedPtr getAlleleMetaData() override { return this->m_allele_meta_data_ptr; }

		virtual uint32_t getForwardCount(std::shared_ptr< Sample > samplePtr, AlleleCountType alleleCountType) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto iter = m_alignment_sample_forward_count_map.find(samplePtr->getName());
			return (iter == m_alignment_sample_forward_count_map.end()) ? 0 : iter->second->getScoreCount(alleleCountType);
		}
		virtual uint32_t getReverseCount(std::shared_ptr< Sample > samplePtr, AlleleCountType alleleCountType) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto iter = m_alignment_sample_reverse_count_map.find(samplePtr->getName());
			return (iter == m_alignment_sample_reverse_count_map.end()) ? 0 : iter->second->getScoreCount(alleleCountType);
		}
		virtual uint32_t getTotalCount(AlleleCountType alleleCountType) override
		{
			uint32_t totalCount = 0;
			for (auto iter : m_alignment_sample_forward_count_map)
			{
				totalCount += iter.second->getTotalCount();
			}
			for (auto iter : m_alignment_sample_reverse_count_map)
			{
				totalCount += iter.second->getTotalCount();
			}
			return totalCount;
		}

		virtual void incrementForwardCount(std::shared_ptr< IAlignment > alignmentPtr, AlleleCountType alleleCountType) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto samplePtr = alignmentPtr->getSample();
			auto iter = m_alignment_sample_forward_count_map.find(samplePtr->getName());
			if (iter == m_alignment_sample_forward_count_map.end())
			{
				/* m_alignment_sample_forward_count_map.emplace(samplePtr->getName(), 1); */
				/* iter = m_alignment_sample_forward_count_map.find(samplePtr->getName()); */
				auto scoreCounterPtr = std::make_shared< ScoreCounter >();
				scoreCounterPtr->incrementScoreCount(alleleCountType);
				m_alignment_sample_forward_count_map.emplace(samplePtr->getName(), scoreCounterPtr);
			}
			else
			{
				iter->second->incrementScoreCount(alleleCountType);
				/* ++iter->second; */
			}
		}

		virtual void incrementReverseCount(std::shared_ptr< IAlignment > alignmentPtr, AlleleCountType alleleCountType) override
		{
			std::lock_guard< std::mutex > lock(m_alignment_count_mutex);
			auto samplePtr = alignmentPtr->getSample();
			auto iter = m_alignment_sample_reverse_count_map.find(samplePtr->getName());
			if (iter == m_alignment_sample_reverse_count_map.end())
			{
				auto scoreCounterPtr = std::make_shared< ScoreCounter >();
				scoreCounterPtr->incrementScoreCount(alleleCountType);
				m_alignment_sample_reverse_count_map.emplace(samplePtr->getName(), std::make_shared< ScoreCounter >());
				/* iter = m_alignment_sample_reverse_count_map.find(samplePtr->getName()); */
			}
			else
			{
				iter->second->incrementScoreCount(alleleCountType);
			}
		}

		virtual void incrementCount(std::shared_ptr< IAlignment > alignmentPtr, AlleleCountType alleleCountType) override
		{
			if (alignmentPtr->isReverseStrand())
			{
				incrementReverseCount(alignmentPtr, alleleCountType);
			}
			else
			{
				incrementForwardCount(alignmentPtr, alleleCountType);
			}
		}

		uint32_t getCommonPrefixSize(IAllele::SharedPtr allelePtr) override
		{
			uint32_t commonIndex = 0;
			uint32_t minSize = (this->getLength() < allelePtr->getLength()) ? this->getLength() : allelePtr->getLength();
			for (uint32_t i = 0; i < minSize; ++i) { if (this->getSequence()[i] != allelePtr->getSequence()[i]) { return i; } }
			return minSize;
		}

		uint32_t getCommonSuffixSize(IAllele::SharedPtr allelePtr) override
		{
			uint32_t commonIndex = 0;
			uint32_t minSize = (this->getLength() < allelePtr->getLength()) ? this->getLength() : allelePtr->getLength();
			for (uint32_t i = 0; i < minSize; ++i) { if (this->getSequence()[(this->getLength() - 1) - i] != allelePtr->getSequence()[(allelePtr->getLength() - 1) - i]) { return i; } }
			return minSize;
		}

	protected:
		Allele() {}

		class ScoreCounter
		{
		public:
			typedef std::shared_ptr< ScoreCounter > SharedPtr;
			ScoreCounter()
			{
				m_score_map[AlleleCountType::NONE] = 0;
				m_score_map[AlleleCountType::AmbiguousPercent] = 0;
				m_score_map[AlleleCountType::SeventyPercent] = 0;
				m_score_map[AlleleCountType::EightyPercent] = 0;
				m_score_map[AlleleCountType::NinteyPercent] = 0;
				m_score_map[AlleleCountType::NinteyFivePercent] = 0;
			}

			void incrementScoreCount(AlleleCountType alleleCountType)
			{
				++m_score_map[alleleCountType];
			}

			uint32_t getScoreCount(AlleleCountType alleleCountType)
			{
				return m_score_map[alleleCountType];
			}

			uint32_t getTotalCount()
			{
				uint32_t count = 0;
				for (auto iter: m_score_map)
				{
					count += iter.second;
				}
				return count;
			}

		private:
			std::unordered_map< AlleleCountType, uint32_t, AlleleCountTypeHash > m_score_map;
		};

		Sequence::SharedPtr m_sequence_ptr;
		AlleleMetaData::SharedPtr m_allele_meta_data_ptr;
		std::mutex m_alignment_count_mutex;
		std::unordered_map< std::string, ScoreCounter::SharedPtr > m_alignment_sample_forward_count_map;
		std::unordered_map< std::string, ScoreCounter::SharedPtr > m_alignment_sample_reverse_count_map;

	};
}

#endif //GRAPHITE_ALLELE_H
