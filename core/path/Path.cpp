#include "Path.h"
#include <functional>

namespace graphite
{
	Path::Path()
	{
	}

	Path::~Path()
	{
	}

	std::vector< IAllele::SharedPtr > Path::getAllelePath()
	{
		return this->m_allele_ptrs;
	}

	uint32_t Path::getPathSWPercent()
	{
		return this->m_sw_percentage;
	}

	IAlignment::SharedPtr Path::getAlignment()
	{
		return this->m_alignment_ptr;
	}

	size_t Path::getHash()
	{
		if (this->m_hash == 0)
		{
			computeAndSetHash();
		}
		return this->m_hash;
	}

	void Path::addAlleleToPath(IAllele::SharedPtr allelePtr)
	{
		this->m_allele_ptrs.emplace_back(allelePtr);
	}

	void Path::setPathSWPercent(uint32_t swPercent)
	{
		this->m_sw_percentage = swPercent;
	}

	void Path::setAlignment(IAlignment::SharedPtr alignmentPtr)
	{
		this->m_alignment_ptr = alignmentPtr;
	}

	void Path::computeAndSetHash()
	{
		// this->m_hash = boost::hash_range(m_allele_ptrs.begin(), m_allele_ptrs.end());
		// std::hash< std::vector< IAllele::SharedPtr > >
		size_t seed = 0;
		for (auto allelePtr : this->m_allele_ptrs)
		{
			seed ^= m_hasher(allelePtr.get()) + 0x9e3779b9 + (seed<<6) + (seed>>2);
			// hash_combine(seed, allelePtr)
		}
		m_hash = seed;
	}

	void Path::setGSSWGraphMapping(std::shared_ptr< gssw_graph_mapping > graphMappingPtr)
	{
		m_graph_mapping_ptr = graphMappingPtr;
	}

	void Path::printLongFormat()
	{
		position startPosition = 0;
		size_t startSoftClipLength = 0;
		size_t endSoftClipLength = 0;
		gssw_node_cigar* nc = this->m_graph_mapping_ptr->cigar.elements;
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
		for (int i = 0; i < this->m_graph_mapping_ptr->cigar.length; ++i, ++nc)
		{
			for (int j = 0; j < nc->cigar->length; ++j)
			{
				if (nc->cigar->elements[j].type == 'S') // softclipping, set softclip variables
				{
					if (i == 0 && j == 0)
					{
						startSoftClipLength = nc->cigar->elements[j].length;
					}
					else
					{
						endSoftClipLength = nc->cigar->elements[j].length;
					}
				}
				cigarString += std::to_string(nc->cigar->elements[j].length) + nc->cigar->elements[j].type;
			}
			if (startPosition == 0)
			{
				startPosition = nc->node->position;
				// referenceStartPosition = startPosition - this->m_reference_ptr->getRegion()->getStartPosition();
			}
			// auto nodeVariantType = static_cast< GenotyperAllele::Type >((long)nc->node->data);
			bool isReferenceNode = (nc->node->id % 2 == 0);
			std::string nodeTypeString = (isReferenceNode) ? "R" : "V";

			cigarString +=  separator;

			// referenceOffsets += std::to_string(referenceStartPosition + refOffset) + separator;
			// size_t totalRefOffset = (referenceStartPosition + refOffset);
			// if (totalRefOffset < this->m_reference_ptr->getRegion()->getEndPosition() - this->m_reference_ptr->getRegion()->getStartPosition())
			// {
				// referenceString += std::string(this->m_reference_ptr->getSequence() + (referenceStartPosition + refOffset), nc->node->ref_len);
			// }

			referenceString += nc->node->ref_seq;
			tracebackString += nc->node->seq;
			nodeSeparatorPositions.emplace_back(nc->node->position + nodeSeparatorOffset);

			if (nc->node->ref_len > nc->node->len) { nodeSeparatorOffset = nc->node->ref_len - nc->node->len; tracebackString += std::string(nc->node->ref_len - nc->node->len, ' '); }
			else if (nc->node->ref_len < nc->node->len) { nodeSeparatorOffset = nc->node->len - nc->node->ref_len; referenceString += std::string(nc->node->len - nc->node->ref_len, ' '); }
			if (i < this->m_graph_mapping_ptr->cigar.length && i != this->m_graph_mapping_ptr->cigar.length - 1)
			{
				referenceString += separator;
				tracebackString += separator;
			}

			// nodeTracebackString += std::string(GenotyperAllele::TypeToString(nodeVariantType)) + separator;
			nodeTracebackString += nodeTypeString + separator;
			nodeOffset += nc->node->len - 1;
			refOffset += nc->node->ref_len;
		}
		nodeTracebackString = (nodeTracebackString.size() > 2) ? nodeTracebackString.substr(0, nodeTracebackString.size() - 2) : nodeTracebackString;
		cigarString = (cigarString.size() > 2) ? cigarString.substr(0, cigarString.size() - 2) : cigarString;

		alignmentString = std::string(this->m_graph_mapping_ptr->position, ' ') + alignmentString.substr(startSoftClipLength);
		// referenceString = std::string(this->m_reference_ptr->getSequence() +(startPosition - this->m_reference_ptr->getRegion()->getStartPosition()), tracebackString.size());

		size_t sepPos = 0;
		while ((sepPos = tracebackString.find(separator, sepPos + 1)) != std::string::npos && sepPos < alignmentString.size())
		{
			alignmentString.insert(sepPos, separator);
		}

		std::string isMapped = (this->m_alignment_ptr->isMapped()) ? "Mapped" : "Unmapped";
		std::string mate = (this->m_alignment_ptr->isFirstMate()) ? "First" : "Second";
		std::string reverseStrand = (this->m_alignment_ptr->isReverseStrand()) ? "Reverse" : "Normal";

		std::string reportString = "";
		std::string eol = "\r\n";
		reportString += "---------------------------------------------------------" + eol;
		reportString += "Score:                " + std::to_string(this->m_graph_mapping_ptr->score) + "/202" + eol;
		reportString += "Node Cigar String:    " + cigarString + eol;
		reportString += "Node Type Traceback:  " + nodeTracebackString + eol;
		reportString += "Reference:            " + referenceString + eol;
		reportString += "Node Traceback:       " + tracebackString + eol;
		reportString += "Alignment:            " + alignmentString + eol;
		reportString += "Alignment(RAW):       " + std::string(this->m_alignment_ptr->getSequence(), this->m_alignment_ptr->getLength()) + eol;
		reportString += "Graph Position:       " + std::to_string(startPosition) + eol;
		reportString += "Mapped State:         " + isMapped + eol;
		reportString += "Mate Order:           " + mate + eol;
		reportString += "Strand:               " + reverseStrand + eol;
		reportString += "Original Map Quality: " + std::to_string(this->m_alignment_ptr->getOriginalMapQuality()) + eol;
		reportString += "---------------------------------------------------------" + eol;

		std::cout << reportString;
	}
}
