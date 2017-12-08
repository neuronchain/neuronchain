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
#pragma once
#include <graphene/chain/protocol/transaction.hpp>

namespace graphene { namespace chain {

   struct packet
   {
       vector<signed_transaction> transactions;
   };
   
   struct signed_packet : public packet
   {
       digest_type                digest()const;
       packet_id_type             id()const;
       void                       sign( const fc::ecc::private_key& signer );
       bool                       validate_signee( const fc::ecc::public_key& expected_signee )const;
       fc::ecc::public_key        signee()const;

       witness_id_type              witness;
       signature_type               witness_signature;
   };

} } // graphene::chain

FC_REFLECT( graphene::chain::packet, (transactions) )
FC_REFLECT_DERIVED( graphene::chain::signed_packet, (graphene::chain::packet), (witness)(witness_signature) )
