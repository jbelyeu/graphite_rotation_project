#include "GSSWMapping.h"

#include <algorithm>
#include <iostream>
#include <thread>

#include <fstream>
#include <sstream>

namespace graphite
{
namespace adjudicator
{
    GSSWMapping::GSSWMapping(std::shared_ptr< gssw_graph_mapping > gsswMappingPtr, IAlignment::SharedPtr alignmentPtr) :
		m_gssw_mapping_ptr(gsswMappingPtr),
		m_alignment_ptr(alignmentPtr),
		m_position(0),
		m_mapped(false)
	{
		uint32_t offset = m_gssw_mapping_ptr->position;
		gssw_node_cigar* nc = m_gssw_mapping_ptr->cigar.elements;
		for (int i = 0; i < m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			if (i == 0)
			{
				this->m_position = nc->node->position + m_gssw_mapping_ptr->position;
			}
			auto allelePtr = ((IAllele*)nc->node->data)->getSharedPtr();
			m_allele_ptrs.push_back(allelePtr);
			m_allele_gssw_nodes_map[allelePtr.get()] = nc->node;
		}
	}

	GSSWMapping::~GSSWMapping()
	{
	}

	std::vector< MappingAlignmentInfo::SharedPtr > GSSWMapping::getMappingAlignmentInfoPtrs(IAdjudicator::SharedPtr adjudicatorPtr)
	{
    static std::mutex smutex;
		std::lock_guard< std::mutex > guard(smutex);
    static std::atomic<bool> fasta_written(false);
    std::string fasta_seq = "";
    uint32_t start_pos = -1;
    bool switchcase = true;

		std::vector< MappingAlignmentInfo::SharedPtr > mappingAlignmentInfoPtrs;
		gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;
		for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			int32_t score = 0;
			uint32_t length = 0;
			uint32_t prefixMatch = 0;
			uint32_t suffixMatch = 0;
			auto allelePtr = ((IAllele*)nc->node->data)->getSharedPtr();
      fasta_seq += std::string(nc->node->seq);
      if (!fasta_written)
      {
        std::string seq = std::string(nc->node->seq);
        std::cout << seq.size();
      }
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				switch (nc->cigar->elements[j].type)
				{
				case 'M':
					score += (adjudicatorPtr->getMatchValue() * nc->cigar->elements[j].length);
					break;
				case 'S':
				case 'X':
					score -= (adjudicatorPtr->getMisMatchValue() * nc->cigar->elements[j].length);
					break;
				case 'I': // I and D are treated the same
				case 'D':
					score -= adjudicatorPtr->getGapOpenValue();
					score -= (adjudicatorPtr->getGapExtensionValue() * (nc->cigar->elements[j].length -1));
					break;
				default:
					break;
				}
				length += nc->cigar->elements[j].length;
			}
			prefixMatch = length;
			suffixMatch = length;
			score = (score < 0) ? 0 : score; // the floor of the mapping score is 0
			auto mappingAlignmentInfo = std::make_shared< MappingAlignmentInfo >(allelePtr, score, length, prefixMatch, suffixMatch);
			mappingAlignmentInfoPtrs.emplace_back(mappingAlignmentInfo);
		}
    if (!fasta_written)
    {
      std::ofstream fasta_file;
      std::stringstream fasta_name;
      //TODO hardcoded output path. Needs fixing.
      fasta_name << "/uufs/chpc.utah.edu/common/home/u1072557/marthlab/" << this->m_gssw_mapping_ptr->graph_start_position <<".fa";
      fasta_file.open(fasta_name.str());
      std::string referenceString;

      //TODO this is a hacky fix for the reference offset. A much better option would be decrementing the read positions,
      //but some (especially the mate positions) become negative because they're before the official graph_start_position
      // for (uint32_t i = 0; i < 498; ++i)
      // {
      //   referenceString += "N";
      // }
      referenceString += fasta_seq;
      fasta_file <<">chr20\n";
      for (uint32_t i = 0; i < referenceString.size(); ++i)
      {
        //TODO I hardcoded the length of the lines for the fasta.
        //We should at the least move this to a class member, or maybe a command-line arg if we're fancy
        if (i > 0 && i%80 == 0)
        {
          fasta_file <<std::endl;
        }
        fasta_file << referenceString[i];
      }

      fasta_written = true;
    }
    // call function to write alignment to sam
    this->outputSAMEntry(this->m_alignment_ptr);



		return mappingAlignmentInfoPtrs;
	}

	/*
	MappingAlignmentInfo::SharedPtr GSSWMapping::getMappingAlignmentInfo(IAllele::SharedPtr allelePtr, std::shared_ptr< IAdjudicator > adjudicatorPtr)
	{
		static std::mutex smutex;
		std::lock_guard< std::mutex > guard(smutex);
		// see above message
		auto iter = m_allele_gssw_nodes_map.find(allelePtr.get());
		if (iter == m_allele_gssw_nodes_map.end() || iter->second->cigar == NULL) { return nullptr; }
		auto gsswNodeCigar = iter->second->cigar;
		int32_t score = 0;
		uint32_t length = 0;
		uint32_t prefixMatch = 0;
		uint32_t suffixMatch = 0;

		for (int j = 0; j < gsswNodeCigar->length; ++j)
		{
			switch (gsswNodeCigar->elements[j].type)
			{
			case 'M':
				if (j == 0) { prefixMatch = gsswNodeCigar->elements[j].length; }
				if (j == gsswNodeCigar->length - 1) { suffixMatch = gsswNodeCigar->elements[j].length; }
				score += (adjudicatorPtr->getMatchValue() * gsswNodeCigar->elements[j].length);
				break;
			case 'X':
				score -= (adjudicatorPtr->getMisMatchValue() * gsswNodeCigar->elements[j].length);
				break;
			case 'I': // I and D are treated the same
			case 'D':
				score -= adjudicatorPtr->getGapOpenValue();
				score -= (adjudicatorPtr->getGapExtensionValue() * (gsswNodeCigar->elements[j].length -1));
				break;
			default:
				break;
			}
			length += gsswNodeCigar->elements[j].length;
		}
		score = (score < 0) ? 0 : score; // the floor of the mapping score is 0
		auto mappingAlignmentInfo = std::make_shared< MappingAlignmentInfo >(allelePtr, score, length, prefixMatch, suffixMatch);
		return  mappingAlignmentInfo;
	}
	*/

	/*
	MappingAlignment::SharedPtr GSSWMapping::getGSSWAlignmentPtrFromAllelePtr(IAllele::SharedPtr allelePtr)
	{
		auto iter = m_allele_mappingalignment_map.find(allelePtr);
		if (iter != m_allele_mappingalignment_map.end())
		{
			return iter->second;
		}
		else
		{
			return nullptr;
		}
	}
	*/

	int GSSWMapping::getMappingScore()
	{
		return this->m_gssw_mapping_ptr->score;
	}

	IAlignment::SharedPtr GSSWMapping::getAlignmentPtr()
	{
		return this->m_alignment_ptr;
	}

	std::vector< IAllele::SharedPtr > GSSWMapping::getAllelePtrs()
	{
		return this->m_allele_ptrs;
	}

	void GSSWMapping::setMapped(bool mapped)
	{
		if (!mapped)
		{
			this->m_allele_incrementor_callback_list.clear();
		}
		this->m_mapped = mapped;
	}

	void GSSWMapping::incrementAlleleCounts()
	{
		if (m_mapped)
		{
			for (auto incrementFunct : this->m_allele_incrementor_callback_list)
			{
				incrementFunct();
			}
		}
	}



  void GSSWMapping::outputSAMEntry(IAlignment::SharedPtr alignmentPtr)
  {
    //TODO: seems odd for the -2 to work. That's probably a bug and may be related to the offset reads
    uint32_t read_offset = this->m_gssw_mapping_ptr->graph_start_position-2;
    int32_t position = this->m_alignment_ptr->getPosition()-this->m_gssw_mapping_ptr->graph_start_position;

    //TODO need a better solution to this problem. I think we should consider clipping
    //off the reads that start before the graph starts, but the cigar string will need
    //to be modified to reflect that. The reads that are being removed here are generally
    //poor matches and may show a problem of some kind in the graph creation.
    if (position <= 0)
    {
      return;
    }

    std::ofstream outfile;
    std::stringstream name;
    //TODO here the filepath is also hardcoded and needs fixing
		name << "/uufs/chpc.utah.edu/common/home/u1072557/marthlab/reads.sam";
		outfile.open(name.str(), std::ios::out | std::ios::app);

    std::vector<std::pair<int32_t, char>> cigar_strings = this->extractCigar();
    std::string stringified_cigar = stringifySimpleCigar(cigar_strings);
    outfile << this->m_alignment_ptr->getID() << "\t";
    outfile << this->m_alignment_ptr->getFlag() << "\t";
    outfile << "chr" << this->m_alignment_ptr->getChrID() +1 << "\t";
    outfile << (this->m_alignment_ptr->getPosition()-read_offset) << "\t";
    outfile << this->m_alignment_ptr->getOriginalMapQuality() << "\t";
    outfile << stringified_cigar << "\t";
    outfile << "chr" << this->m_alignment_ptr->getRefNext() << "\t";
    int32_t pos_next = (this->m_alignment_ptr->getPositionNext()-read_offset);

    //TODO this is just to avoid the samtools error for mates mapped to 0. Probably shouldn't be done.
    if (pos_next <= 0)
    {
      pos_next = 1;
    }
    outfile << pos_next << "\t";
    outfile << this->m_alignment_ptr->getTemplateLength() << "\t";
    outfile << alignmentPtr->getSequence() << "\t";
    //TODO: this is hardcoded QUAL info. I couldn't find a way to easily parse it out for the SAM. so I hard-coded it to speed things up.
    outfile << "000770<<B<B7B<<7<BBB<B0<B<0<<00BB0<<BBBB<BBBBB0<B<<BB<<<<B0BB0<<BB<0BBBB<BB<BB00<<BBBBB00<007<B<BB00BB<BB<BBBBB<00<<<B<BBB07<B	AS:i:126	MC:Z:126M	MQ:i:0	XA:Z:chrX,-156022948,126M,0;chr1,+187516,126M,0;chr16,+16677,126M,0;chr12,+17109,126M,0;chr15,-101973841,126M,0;chr2,-113596317,126M,1;chr9,+17105,126M,2;chr12_GL877875v1_alt,+7109,126M,0;	XS:i:126	MD:Z:126	NM:i:0RG:Z:ERR894730\n";
  }

  std::string GSSWMapping::stringifySimpleCigar(std::vector<std::pair<int32_t, char>> simple_cigar)
  {
    std::stringstream ss;
    for (auto iter = simple_cigar.begin(); iter != simple_cigar.end(); ++iter)
    {
      ss << iter->first << iter->second;
    }
    return ss.str();
  }

  std::vector<std::pair<int32_t, char>> GSSWMapping::extractCigar()
  {
    std::vector<std::pair<int32_t, char>> simple_cigar;
    for (uint32_t i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i)
    {
      for (uint32_t j = 0; j < this->m_gssw_mapping_ptr->cigar.elements[i].cigar->length; ++j)
      {
        simple_cigar.emplace_back(this->m_gssw_mapping_ptr->cigar.elements[i].cigar->elements[j].length,
          this->m_gssw_mapping_ptr->cigar.elements[i].cigar->elements[j].type);
      }
    }
    return simple_cigar;
  }

	void GSSWMapping::addAlleleCountCallback(std::function< void () > functor)
	{
		this->m_allele_incrementor_callback_list.emplace_back(functor);
	}

	void GSSWMapping::printSimpleMapping()
	{
		gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;
		std::string nodeTypesString = "";
		std::string positions = "";
		std::string alignmentString = std::string(this->m_alignment_ptr->getSequence(), this->m_alignment_ptr->getLength());
		std::string referenceString = "";
		std::string tracebackString = "";
		std::string separator = "||";
		for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			std::string nodeAlignmentSequence = std::string(nc->node->seq, nc->node->len);
			std::string nodeReferenceSequence = std::string(nc->node->ref_seq, nc->node->ref_len);
			positions += std::to_string(nc->node->position) + " ";
			size_t tracebackOffset = 0;
			int refSizeDiff = nc->node->ref_len - nc->node->len;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				uint32_t cigLen = nc->cigar->elements[j].length;
				std::string cigarTracebackString = std::string(nc->node->seq + tracebackOffset, cigLen);
				switch (nc->cigar->elements[j].type)
				{
				case 'S':
				case 'X':
					std::transform(cigarTracebackString.begin(), cigarTracebackString.end(), cigarTracebackString.begin(), ::tolower);
					break;
				case 'M':
				case 'I':
				case 'D':
				default:
					break;
				}
				tracebackString += cigarTracebackString;
			}
			tracebackString += separator;
			nodeTypesString += (nc->node->id % 2 == 0) ? "R" + separator : "V" + separator;
		}
		std::string currentMapping = (this->m_mapped) ? "Mapped" : "Unmapped";
		std::cout << "----------------------------" << std::endl;
		std::cout << "Node Types: " << nodeTypesString.substr(0, nodeTypesString.size() - 2) << std::endl;
		std::cout << "Score:                " << std::to_string(this->m_gssw_mapping_ptr->score) << std::endl;
		std::cout << "Current Mapped State: " << currentMapping << std::endl;
		std::cout << "Positions: " << std::endl;
		std::cout << "Reference: " << referenceString.substr(0, referenceString.size() - 2) << std::endl;
		std::cout << "Traceback: " << tracebackString.substr(0, tracebackString.size() - 2) << std::endl;
		std::cout << "Alignment: " << alignmentString.substr(0, alignmentString.size() - 2) << std::endl;
		std::cout << "----------------------------" << std::endl;
	}

	void GSSWMapping::printMapping()
	{
		position startPosition = 0;
		size_t startSoftClipLength = 0;
		size_t endSoftClipLength = 0;
		gssw_node_cigar* nc = this->m_gssw_mapping_ptr->cigar.elements;
		std::string separator = "||";
		std::string nodeTracebackString = "";
		std::string tracebackString = "";
		std::string cigarString = "";
		std::string alignmentString = this->m_alignment_ptr->getSequence();

		std::vector< position > nodeSeparatorPositions;
		// position referenceStartPosition = 0;
		std::string referenceString = "";
		size_t nodeOffset = 0;
		size_t refOffset = 0;
		size_t nodeSeparatorOffset = 0;

		std::string referenceOffsets;
		std::string tmpAlignmentString;
		std::string tmpTracebackString;
		std::string tmpReferenceString;
		size_t alignmentOffset = 0;
		for (int i = 0; i < this->m_gssw_mapping_ptr->cigar.length; ++i, ++nc)
		{
			size_t offset = 0;
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				switch (nc->cigar->elements[j].type)
				{
				case 'S':
					if (i == 0 && j == 0)
					{
						startSoftClipLength = nc->cigar->elements[j].length;
						alignmentOffset += startSoftClipLength;
						// tmpAlignmentString += std::string(nc->cigar->elements[j].length, ' ');
					}
					else
					{
						endSoftClipLength = nc->cigar->elements[j].length;
						alignmentString.erase(alignmentString.size() - endSoftClipLength - 1, endSoftClipLength);
					}
					// tmpReferenceString += std::string(nc->node->ref_seq + offset, nc->cigar->elements[j].length);
					break;
				case 'M':
				case 'X':
					tmpAlignmentString += std::string(alignmentString.c_str() + alignmentOffset, nc->cigar->elements[j].length);
					tmpTracebackString += std::string(nc->node->seq + offset, nc->cigar->elements[j].length);
					// tmpReferenceString += std::string(nc->node->ref_seq + offset, nc->cigar->elements[j].length);
					offset += nc->cigar->elements[j].length;
					alignmentOffset += nc->cigar->elements[j].length;
					break;
				case 'I':
					tmpAlignmentString += std::string(alignmentString.c_str() + alignmentOffset, nc->cigar->elements[j].length);
					tmpTracebackString += std::string(nc->cigar->elements[j].length, '-');
					// tmpReferenceString += std::string(nc->node->ref_seq + offset, nc->cigar->elements[j].length);
					offset += nc->cigar->elements[j].length;
					alignmentOffset += nc->cigar->elements[j].length;
					break;
				case 'D':
					tmpAlignmentString += std::string(alignmentString.c_str() + alignmentOffset, nc->cigar->elements[j].length);
					// tmpReferenceString += std::string(nc->cigar->elements[j].length, '-');
					tmpTracebackString += std::string(nc->node->seq + offset, nc->cigar->elements[j].length);
					offset += nc->cigar->elements[j].length;
					alignmentOffset += nc->cigar->elements[j].length;
					break;
				}
				cigarString += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
			}
			if (startPosition == 0)
			{
				startPosition = nc->node->position;
			}
			std::string nodeTypeString = (nc->node->id % 2 == 0) ? "R" : "V";

			cigarString +=  separator;
			tmpAlignmentString += separator;
			tmpTracebackString += separator;
			tmpReferenceString += std::string(nc->node->ref_seq, nc->node->ref_len) + separator;

			referenceString += nc->node->ref_seq;
			tracebackString += nc->node->seq;
			nodeSeparatorPositions.emplace_back(nc->node->position + nodeSeparatorOffset);

			if (nc->node->ref_len > nc->node->len) { nodeSeparatorOffset = nc->node->ref_len - nc->node->len; tracebackString += std::string(nc->node->ref_len - nc->node->len, ' '); }
			else if (nc->node->ref_len < nc->node->len) { nodeSeparatorOffset = nc->node->len - nc->node->ref_len; referenceString += std::string(nc->node->len - nc->node->ref_len, ' '); }
			if (i < this->m_gssw_mapping_ptr->cigar.length && i != this->m_gssw_mapping_ptr->cigar.length - 1)
			{
				referenceString += separator;
				tracebackString += separator;
			}

			nodeTracebackString += nodeTypeString + separator;
			nodeOffset += nc->node->len - 1;
			refOffset += nc->node->ref_len;
		}
		nodeTracebackString = (nodeTracebackString.size() > 2) ? nodeTracebackString.substr(0, nodeTracebackString.size() - 2) : nodeTracebackString;
		cigarString = (cigarString.size() > 2) ? cigarString.substr(0, cigarString.size() - 2) : cigarString;
		tmpAlignmentString = (tmpAlignmentString.size() > 2) ? tmpAlignmentString.substr(0, tmpAlignmentString.size() - 2) : tmpAlignmentString;
		tmpTracebackString = (tmpTracebackString.size() > 2) ? tmpTracebackString.substr(0, tmpTracebackString.size() - 2) : tmpTracebackString;
		tmpReferenceString = (tmpReferenceString.size() > 2) ? tmpReferenceString.substr(0, tmpReferenceString.size() - 2) : tmpReferenceString;

		alignmentString = std::string(this->m_gssw_mapping_ptr->position, ' ') + alignmentString.substr(startSoftClipLength);

		size_t sepPos = 0;
		while ((sepPos = tracebackString.find(separator, sepPos + 1)) != std::string::npos && sepPos < alignmentString.size())
		{
			alignmentString.insert(sepPos, separator);
		}

		std::string isMapped = (this->m_alignment_ptr->isMapped()) ? "Mapped" : "Unmapped";
		std::string mate = (this->m_alignment_ptr->isFirstMate()) ? "First" : "Second";
		std::string reverseStrand = (this->m_alignment_ptr->isReverseStrand()) ? "Reverse" : "Normal";
		std::string currentMapping = (this->m_mapped) ? "Mapped" : "Unmapped";

		std::string reportString = "";
		std::string eol = "\r\n";
		reportString += "---------------------------------------------------------" + eol;
		reportString += "Score:                " + std::to_string(this->m_gssw_mapping_ptr->score) + eol;
		reportString += "Node Cigar String:    " + cigarString + eol;
		reportString += "Node Type Traceback:  " + nodeTracebackString + eol;
		reportString += "Reference:            " + referenceString + eol;
		// reportString += "Reference (new):      " + tmpReferenceString + eol;
		// reportString += "Alignment (new):      " + tmpAlignmentString + eol;
		// reportString += "Traceback (new):      " + tmpTracebackString + eol;
		reportString += "Node Traceback:       " + tracebackString + eol;
		reportString += "Alignment:            " + alignmentString + eol;
		reportString += "Alignment(RAW):       " + std::string(this->m_alignment_ptr->getSequence(), this->m_alignment_ptr->getLength()) + eol;
		reportString += "Alignment Name:       " + this->m_alignment_ptr->getID() + eol;
		reportString += "Alignment Position:   " + std::to_string(this->m_alignment_ptr->getPosition()) + eol;
		reportString += "Graph Position:       " + std::to_string(startPosition) + eol;
		reportString += "Mapping ID:           " + std::to_string(this->m_id) + eol;
		reportString += "Current Mapped State: " + currentMapping + eol;
		reportString += "Prev. Mapped State:   " + isMapped + eol;
		reportString += "Mate Order:           " + mate + eol;
		reportString += "Strand:               " + reverseStrand + eol;
		reportString += "Original Map Quality: " + std::to_string(this->m_alignment_ptr->getOriginalMapQuality()) + eol;
		reportString += "---------------------------------------------------------" + eol;

		std::cout << reportString;
	}

}
}
