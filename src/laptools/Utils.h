/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/* MIT License                                                                     */
/*                                                                                 */
/* Copyright (c) 2022 Nariaki Tateiwa <n-tateiwa@kyudai.jp>                        */
/*                                                                                 */
/* Permission is hereby granted, free of charge, to any person obtaining a copy    */
/* of this software and associated documentation files (the "Software"), to deal   */
/* in the Software without restriction, including without limitation the rights    */
/* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       */
/* copies of the Software, and to permit persons to whom the Software is           */
/* furnished to do so, subject to the following conditions:                        */
/*                                                                                 */
/* The above copyright notice and this permission notice shall be included in all  */
/* copies or substantial portions of the Software.                                 */
/*                                                                                 */
/* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      */
/* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        */
/* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     */
/* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          */
/* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   */
/* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   */
/* SOFTWARE.                                                                       */
/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

/**@file    Utils.h
 * @brief   Utilities for laptools.
 * @author  Nariaki Tateiwa
 *
 *
 *
 */

/*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/


#ifndef __LAPTOOLS_UTILS_H__
#define __LAPTOOLS_UTILS_H__

#include <vector>
#include <unordered_map>

namespace LapTools
{

template <typename T>
class CountArray
{

public:
   int         length;                 ///< length of array
   int         index;                  ///< index of current basis hash
   std::vector<T> array;               ///<
   std::unordered_map<T, int> counter; ///< key is elment, value is number of elment
   // int         mostFreqElm;            ///< most frequency element

   ///
   /// @brief defalt constructor
   ///
   CountArray() : index(0) {};

   ///
   /// @brief Constructor
   ///
   void resize(
         int inLength
         )
   {
      index = 0;
      length = inLength;
      array.resize(length);
      counter.clear();
      for( T i = length; i > 0; --i )
      {
         array[index++] = i;
         counter.emplace(i, 1);
      }
      index = 0;
   }

   ///
   /// @breif clear data
   ///
   void clear(
         )
   {
      resize(length);
   }

   ///
   /// @brief push new element
   /// @param[in] elm element
   /// @details add/delete element into array, and update counter
   ///
   void push(
         T elm
         )
   {
      T deleteElm = array[index];
      array[index++] = elm;
      if( index == length ) index = 0;
      if( counter[deleteElm] > 1) counter[deleteElm]--;
      else counter.erase(deleteElm);
      if( counter.count(elm) > 0 ) counter[elm]++;
      else counter.emplace(elm, 1);
   }

   ///
   /// @brief get number of most frequency elment in array
   ///
   T mostFrequentCount(
         )
   {
      int __mostFrequentCount = -1;
      for( auto pair : counter )
      {
         if( pair.second > __mostFrequentCount )
            __mostFrequentCount = pair.second;
      }
      return __mostFrequentCount;
   }

}; // class CountArray


}  // namespace LapTools


#endif // __LAPTOOLS__UTILS_H__
