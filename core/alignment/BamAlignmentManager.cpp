#include "BamAlignmentManager.h"
#include "BamAlignmentReader.h"
#include "SamtoolsAlignmentReader.h"
#include "AlignmentList.h"
#include "core/util/ThreadPool.hpp"

#include <functional>
#include <boost/filesystem.hpp>
#include <unordered_set>
#include <algorithm>
#include <deque>
#include <cstdlib>

namespace graphite
{
	BamAlignmentManager::BamAlignmentManager(const std::vector< Sample::SharedPtr >& samplePtrs, Region::SharedPtr regionPtr, bool excludeDuplicateReads) :
		m_region_ptr(regionPtr),
		m_sample_ptrs(samplePtrs),
		m_loaded(false),
        m_exclude_duplicate_reads(excludeDuplicateReads)
	{
    }

	BamAlignmentManager::~BamAlignmentManager()
	{
	}

	IAlignmentList::SharedPtr BamAlignmentManager::getAlignmentsInRegion(Region::SharedPtr regionPtr)
	{
		position startPosition = regionPtr->getStartPosition();
		position endPosition = regionPtr->getEndPosition();

		auto lowerBound = std::lower_bound(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(), nullptr, [startPosition](const IAlignment::SharedPtr& alignmentPtr, const IAlignment::SharedPtr& ignore) {
				return startPosition > alignmentPtr->getPosition() + alignmentPtr->getLength();
			});
		auto upperBound = std::upper_bound(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(), nullptr, [endPosition](const IAlignment::SharedPtr& ignore, const IAlignment::SharedPtr& alignmentPtr) {
				return alignmentPtr->getPosition() > endPosition;
			});

		std::vector< IAlignment::SharedPtr > alignmentPtrs;
		alignmentPtrs.insert(alignmentPtrs.begin(), lowerBound, upperBound);
		return std::make_shared< AlignmentList >(alignmentPtrs);
	}

	void BamAlignmentManager::asyncLoadAlignments(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		if (this->m_loaded) { return; }
		std::unordered_set< std::string > usedPaths;
		for (auto samplePtr : this->m_sample_ptrs)
		{
			if (usedPaths.find(samplePtr->getPath()) != usedPaths.end()) { continue; }
			usedPaths.emplace(samplePtr->getPath());
			this->m_loading_thread_ptrs.emplace_back(std::make_shared< std::thread >(&BamAlignmentManager::loadBam, this, samplePtr->getPath(), variantManagerPtr, variantPadding));
		}
	}

	std::vector< Region::SharedPtr > BamAlignmentManager::getRegionsContainingVariantsWithPadding(IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		std::vector< Region::SharedPtr > regionPtrs;

		auto variantListPtr = variantManagerPtr->getVariantsInRegion(this->m_region_ptr);
		std::vector< IVariant::SharedPtr > variantPtrs;
		IVariant::SharedPtr variantPtr;
		while (variantListPtr->getNextVariant(variantPtr))
		{
			variantPtrs.emplace_back(variantPtr);
		}

		Region::SharedPtr regionPtr = nullptr;
		for (auto i = 0; i < variantPtrs.size(); ++i)
		{
			// std::cout << "variant: " << i << std::endl;
			position startPosition = ((variantPtrs[i]->getPosition() - (variantPadding * 2)) > 0) ? (variantPtrs[i]->getPosition() - (variantPadding * 2)) : 0;
			position endPosition = (variantPtrs[i]->getPosition() + variantPtrs[i]->getMaxAlleleSize() + (variantPadding * 2));
			auto j = i + 1;
			while (j < (variantPtrs.size() - 1) && variantPtrs[j]->getPosition() < endPosition)
			{
				endPosition = (variantPtrs[j]->getPosition() + variantPtrs[j]->getMaxAlleleSize() + (variantPadding * 2));
				++j;
			}
			i = j - 1;
			regionPtr = std::make_shared< Region >(this->m_region_ptr->getReferenceID(), startPosition, endPosition);
			regionPtrs.push_back(regionPtr);
		}
		return regionPtrs;
	}

	void BamAlignmentManager::loadBam(const std::string bamPath, IVariantManager::SharedPtr variantManagerPtr, uint32_t variantPadding)
	{
		ThreadPool::Instance()->start();
		auto regionPtrs = getRegionsContainingVariantsWithPadding(variantManagerPtr, variantPadding);

		std::deque< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;

		// std::vector< std::shared_ptr< std::future< std::vector< IAlignment::SharedPtr > > > > futureFunctions;


		std::vector< SamtoolsAlignmentReader::SharedPtr > bamAlignmentReaders;
		std::unordered_set< std::string > regionStrings;
		for (auto regionPtr : regionPtrs)
		{
			position bamLastPosition = SamtoolsAlignmentReader::GetLastPositionInBam(bamPath, regionPtr);
			// std::cout << "end position: " << regionPtr->getReferenceID() << " " << regionPtr->getStartPosition() <<  " " << regionPtr->getEndPosition() << std::endl;
			// position lastRegionPosition = (bamLastPosition < regionPtr->getEndPosition()) ? bamLastPosition : regionPtr->getEndPosition();
			auto bamAlignmentReaderPtr = std::make_shared< SamtoolsAlignmentReader >(bamPath);
			auto funct = std::bind(&SamtoolsAlignmentReader::loadAlignmentsInRegion, bamAlignmentReaderPtr, regionPtr, this->m_exclude_duplicate_reads);
			auto future = ThreadPool::Instance()->enqueue(funct);
			futureFunctions.push_back(future);
			bamAlignmentReaders.push_back(bamAlignmentReaderPtr);
		}

		std::unordered_set< std::string > alignmentSet;
		while (!futureFunctions.empty())
		{
			auto futureFunct = futureFunctions.front();
			futureFunctions.pop_front();
			if (futureFunct->wait_for(std::chrono::milliseconds(100)) == std::future_status::ready)
			{
				std::lock_guard< std::mutex > lockGuard(m_alignment_ptrs_lock);
				auto loadedAlignmentPtrs = futureFunct->get();
				for (auto& alignment : loadedAlignmentPtrs)
				{
					if (alignmentSet.find(alignment->getID()) == alignmentSet.end())
					{
						alignmentSet.insert(alignment->getID());
						this->m_alignment_ptrs.push_back(alignment);
					}
				}
				// std::cout << "added: " << loadedAlignmentPtrs.size() << std::endl;
				continue;
			}
			else
			{
				futureFunctions.emplace_back(futureFunct);
			}
			// futureFunct->wait();
		}

		/*
		std::unordered_set< std::string > alignmentSet;
		for (auto futureFunct : futureFunctions)
		{
			futureFunct->wait();
		}

		{
			std::lock_guard< std::mutex > lockGuard(m_alignment_ptrs_lock);
			for (auto bamAlignmentReader : bamAlignmentReaders)
			{
				for (auto& alignment : bamAlignmentReader->getBamAlignments())
				{
					if (alignmentSet.find(alignment->getID()) == alignmentSet.end())
					{
						alignmentSet.insert(alignment->getID());
						this->m_alignment_ptrs.push_back(alignment);
					}
				}
			}
		}
		*/
		// std::cout << "finished loading alignments: " << this->m_alignment_ptrs.size() << std::endl;
		this->m_loaded = true;
	}

	std::vector< Sample::SharedPtr > BamAlignmentManager::getSamplePtrs() { return m_sample_ptrs; }

	void BamAlignmentManager::waitForAlignmentsToLoad()
	{
		for (auto threadPtr : this->m_loading_thread_ptrs)
		{
			threadPtr->join();
		}
		// sort the alignments since they could be from different bam files
		std::sort(this->m_alignment_ptrs.begin(), this->m_alignment_ptrs.end(),
				  [] (const IAlignment::SharedPtr& a, const IAlignment::SharedPtr& b)
				  {
					  return a->getPosition() < b->getPosition();
				  });

		for (auto alignmentPtr : this->m_alignment_ptrs)
		{
			std::cout << alignmentPtr->getPosition() << " " << alignmentPtr->getSequence() << std::endl;
		}
	}

	void BamAlignmentManager::releaseResources()
	{
	}

	void BamAlignmentManager::processMappingStatistics()
	{
		for (auto alignmentPtr :  this->m_alignment_ptrs)
		{
			auto alignmentMappingMutexPtr = alignmentPtr->getMappingMutex();
			std::lock_guard< std::recursive_mutex > lock(*alignmentMappingMutexPtr); // make sure the alignmentMapping isn't set during this
			if (auto alignmentMappingWPtr = alignmentPtr->getMapping().lock()) // check if the alignment has already been aligned previously
			{
			}
		}
	}
}
