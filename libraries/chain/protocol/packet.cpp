/*
 * Copyright (c) 2015 Cryptonomex, Inc., and contributors.
 *
 * The MIT License
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
#include <graphene/chain/protocol/packet.hpp>
#include <fc/io/raw.hpp>
#include <fc/bitutil.hpp>
#include <algorithm>

namespace graphene { namespace chain {
   digest_type signed_packet::digest()const
   {    
      std::vector<digest_type> hashes;
      for (const auto& tx : transactions)
        hashes.push_back(tx.digest());

      digest_type::encoder enc;
      fc::raw::pack( enc, hashes );
      return enc.result();
   }

   graphene::chain::packet_id_type signed_packet::id() const
   {
      auto h = digest();
      packet_id_type result;
      memcpy(result._hash, h._hash, std::min(sizeof(result), sizeof(h)));
      return result;
   }

   fc::ecc::public_key signed_packet::signee()const
   {
      return fc::ecc::public_key( witness_signature, signed_packet::digest(), true/*enforce canonical*/ );
   }

   void signed_packet::sign( const fc::ecc::private_key& signer )
   {
      witness_signature = signer.sign_compact( digest() );
   }

   bool signed_packet::validate_signee( const fc::ecc::public_key& expected_signee )const
   {
      return signee() == expected_signee;
   }

   /*checksum_type signed_block::calculate_merkle_root()const
   {
      if( transactions.size() == 0 ) 
         return checksum_type();

      vector<digest_type> ids;
      ids.resize( transactions.size() );
      for( uint32_t i = 0; i < transactions.size(); ++i )
         ids[i] = transactions[i].merkle_digest();

      vector<digest_type>::size_type current_number_of_hashes = ids.size();
      while( current_number_of_hashes > 1 )
      {
         // hash ID's in pairs
         uint32_t i_max = current_number_of_hashes - (current_number_of_hashes&1);
         uint32_t k = 0;

         for( uint32_t i = 0; i < i_max; i += 2 )
            ids[k++] = digest_type::hash( std::make_pair( ids[i], ids[i+1] ) );

         if( current_number_of_hashes&1 )
            ids[k++] = ids[i_max];
         current_number_of_hashes = k;
      }
      return checksum_type::hash( ids[0] );
   }*/

} }
