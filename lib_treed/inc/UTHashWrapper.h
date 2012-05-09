/*
  This file is part of INDDGO.

  Copyright (C) 2012, Oak Ridge National Laboratory 

  This product includes software produced by UT-Battelle, LLC under Contract No. 
  DE-AC05-00OR22725 with the Department of Energy. 

  This program is free software; you can redistribute it and/or modify
  it under the terms of the New BSD 3-clause software license (LICENSE). 
  
  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
  LICENSE for more details.

  For more information please contact the INDDGO developers at: 
  inddgo-info@googlegroups.com

*/

#ifndef __UTHASH_WRAPPER__
#define __UTHASH_WRAPPER__
/*
 * class wrapper for uthash. This wrapper makes easier to add and find elements from uthash.
 */
class UTHashWrapper
{
  private:
	TDSolution *&head;
	int masklen;

  public:
	UTHashWrapper(TDSolution *&head, int masklen) :
		head(head), masklen(masklen)
	{
	};
        
        ~UTHashWrapper()
        {
        };

        inline void hash_add(bigint_t& mask, int value)
        {
            TDSolution *added_set;
            added_set = new TDSolution(masklen);
            *(added_set->mask) = mask;
            added_set->value = value;
            HASH_ADD_KEYPTR(hh, head, added_set->mask->words,
                            masklen * sizeof(BIGINT_WORD), added_set);
        }
        ;

        inline TDSolution* hash_find(bigint_t& mask)
        {
            TDSolution *hash_ptr = NULL;
            HASH_FIND(hh, head, mask.words, masklen * sizeof(BIGINT_WORD), hash_ptr);
            return hash_ptr;
        }
        ;

        int size()
        {
            return HASH_COUNT(head);
        }
        ;

        inline void hash_delete()
        {
            TDSolution *current_set;
            TDSolution *tmp;
            HASH_ITER(hh, head, current_set, tmp)
            {
                HASH_DEL(head,current_set);
                delete current_set;
            }
            head=NULL;
        };
        
        // This method can be used to access the head when iterate
        // through the hash table.
        TDSolution* get_head() const
        {
            return head;
        };
};
#endif  /* __UTHASH_WRAPPER__ */

