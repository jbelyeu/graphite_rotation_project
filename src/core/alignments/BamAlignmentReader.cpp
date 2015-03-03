#include "BamAlignmentReader.h"
#include "BamAlignment.h"

namespace gwiz
{

	BamAlignmentReader::BamAlignmentReader(const std::string& bamPath) :
		m_bam_path(bamPath), m_average_bam_read_length(0)
	{
		if (!m_bam_reader.Open(m_bam_path))
		{
			throw "Unable to open bam file";
		}
		m_bam_reader.GetHeaderText();
		m_bam_reader.GetReferenceData();
		setAverageBamReadLength();
	}

	BamAlignmentReader::~BamAlignmentReader()
	{
		m_bam_reader.Close();
	}

	void BamAlignmentReader::setAverageBamReadLength()
	{
		IAlignment::SharedPtr alignmentPtr;
		if (getNextAlignment(alignmentPtr))
		{
			m_average_bam_read_length = alignmentPtr->getLength();
			this->m_bam_reader.Rewind();
		}
	}

	bool BamAlignmentReader::getNextAlignment(IAlignment::SharedPtr& alignmentPtr)
	{
		BamAlignmentPtr bamAlignmentPtr = std::make_shared< BamTools::BamAlignment >();
		bool returnValue = this->m_bam_reader.GetNextAlignment(*bamAlignmentPtr);
		alignmentPtr = std::make_shared< BamAlignment >(bamAlignmentPtr);

		return returnValue;
	}

	void BamAlignmentReader::setRegion(Region::SharedPtr region)
	{
		int refID = this->m_bam_reader.GetReferenceID(region->getReferenceID());
		this->m_bam_reader.SetRegion(refID, region->getStartPosition(), refID, region->getEndPosition());
	}

}