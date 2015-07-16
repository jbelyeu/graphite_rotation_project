#ifndef GWIZ_EQUIVALENTALLELE_H
#define GWIZ_EQUIVALENTALLELE_H

#include "Allele.h"

namespace gwiz
{
	/*
	 * The purpose of this allele class is to contain alleles
	 * that have the same sequence but different padding and
	 * therefore can be part of multiple VCF files. The reason
	 * that the alleles are contained all within this allele
	 * is so that when it is adjudicated the adjudication can
	 * be propagated to all alleles in m_allele_ptrs.
	 */
	class EquivalentAllele : public Allele
	{
	public:
	    EquivalentAllele(Sequence::SharedPtr sequence) : Allele(sequence) {}
	    EquivalentAllele(const std::string& seqString, std::vector< IAllele::SharedPtr > allelePtrs) : Allele(seqString)
		{
			m_allele_ptrs = allelePtrs;
			/* std::cout << "EquivalentAllele: " << m_sequence_ptr->getLength() << std::endl; */
		}
		~EquivalentAllele() {}

		void addAllele(IAllele::SharedPtr allelePtr) { this->m_allele_ptrs.emplace_back(allelePtr); }
		std::vector< IAllele::SharedPtr > getAllAlleles() { return this->m_allele_ptrs; }
	private:
		std::vector< IAllele::SharedPtr > m_allele_ptrs;
		/* Sequence::SharedPtr m_sequence_ptr; */
	};
}

#endif //GWIZ_EQUIVALENTALLELE_H
