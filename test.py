# /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
# /* MIT License                                                                     */
# /*                                                                                 */
# /* Copyright (c) 2022 Nariaki Tateiwa <n-tateiwa@kyudai.jp>                        */
# /*                                                                                 */
# /* Permission is hereby granted, free of charge, to any person obtaining a copy    */
# /* of this software and associated documentation files (the "Software"), to deal   */
# /* in the Software without restriction, including without limitation the rights    */
# /* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell       */
# /* copies of the Software, and to permit persons to whom the Software is           */
# /* furnished to do so, subject to the following conditions:                        */
# /*                                                                                 */
# /* The above copyright notice and this permission notice shall be included in all  */
# /* copies or substantial portions of the Software.                                 */
# /*                                                                                 */
# /* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR      */
# /* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,        */
# /* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE     */
# /* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER          */
# /* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,   */
# /* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE   */
# /* SOFTWARE.                                                                       */
# /* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

# /**@file    test.py
#  * @brief   test for CMAP-LAP and CMAP-DeepBKZ
#  * @author  Nariaki Tateiwa
#  *
#  *
#  *
#  */

# /*---+----1----+----2----+----3----+----4----+----5----+----6----+----7----+----8----+----9----+----0----+----1----+----2*/



import sys
import time
import argparse
import subprocess
import collections


Test = collections.namedtuple("Test", "name command")


def main(args):

    # generate test command list
    test_list = gen_test_list(args)

    # test
    if args.only_commands:
        num_test = len(test_list)
        for i, test in enumerate(test_list):
            print(f"({i+1}/{num_test} {test.name}) {test.command}")
            print()
    else:
        if args.long:
            timeout = 300 * 3;
        else:
            timeout = 30 * 3;
        test_commands(test_list, args, timeout)


def gen_test_list(args):
    """generate test commands

    Parameters
    ----------
    args : argparse.Namespace

    Returns
    -------
    test_list : list of Test
    """

    #
    # setting
    #
    sth = args.number_of_parallels
    instance = args.instance
    base_setting_files = "./settings/test/debug.set"
    if args.long:
        base_setting_files += " ./settings/test/long.set"
    else:
        base_setting_files += " ./settings/test/short.set"

    test_list = list()

    #
    # test for fcmapdeepbkz
    #
    if not args.no_test_shared:
        binary = "./bin/fcmapdeepbkz"
        appendix_setting_files_list = [
            (
                "thread -- normal",
                ""
            ),
            (
                "thread -- with lower bound of approximated factor",
                "loweralpha.set"
            )
        ]

        for test_name, appendix_setting_files in appendix_setting_files_list:
            num_setting_files = len(base_setting_files.split()) + len(appendix_setting_files.split())
            appendix_setting_files = " ".join([f"./settings/test/{sfile}" for sfile in appendix_setting_files.split()])
            test_list.append(
                Test(
                    test_name,
                    f"{binary} ./settings/default.set {instance} -ntpr 1 -sth {sth} -o {num_setting_files} {base_setting_files} {appendix_setting_files}",
                )
            )

    #
    # test for paracmaptest
    #
    if not args.no_test_distributed:
        binary = "./bin/paracmapdeepbkz"
        appendix_setting_files_list = [
            (
                "mpi -- noamal",
                ""
            ),
            (
                "mpi-- with lower bound of approximated factor",
                "loweralpha.set"
            ),
        ]
        for test_name, appendix_setting_files in appendix_setting_files_list:
            num_setting_files = len(base_setting_files.split()) + len(appendix_setting_files.split())
            appendix_setting_files = " ".join([f"./settings/test/{sfile}" for sfile in appendix_setting_files.split()])
            test_list.append(
                Test(
                    test_name,
                    f"mpirun -n {sth} {binary} ./settings/default.set {instance} -ntpr 1 -o {num_setting_files} {base_setting_files} {appendix_setting_files}"
                )
            )

    return test_list


def test_commands(test_list, args, timeout):
    """execute commands and display results

    Parameters
    ----------
    test_list : list of Test
    args : argparse.Namespace
    timeout: float
    """
    num_test = len(test_list)
    results = list()
    times = list()
    for i, test in enumerate(test_list):
        sys.stdout.write(f"({i+1}/{num_test} {test.name})[{cyan_str('Checking..')}] {test.command}")
        sys.stdout.flush()
        start_time = time.time()
        try:
            if args.verbose:
                stdout, stderr = None, None
            else:
                stdout, stderr = subprocess.DEVNULL, subprocess.DEVNULL
            subprocess.run(test.command, check=True, shell=True, stdout=stdout, stderr=stderr, timeout=timeout)
        except subprocess.TimeoutExpired:
            print(f"\r({i+1}/{num_test} {test.name})[{red_str('Timeout')}] time {time.time()-start_time:.2f}s  {test.command}")
            results.append("Timeout")
            times.append(time.time()-start_time)
            continue
        except subprocess.CalledProcessError:
            print(f"\r({i+1}/{num_test} {test.name})[{red_str('NG')}] time {time.time()-start_time:.2f}s  {test.command}")
            results.append("NG")
            times.append(time.time()-start_time)
            continue
        print(f"\r({i+1}/{num_test} {test.name})[{green_str('OK')}] time {time.time()-start_time:.2f}s  {test.command}")
        results.append("OK")
        times.append(time.time()-start_time)
    if args.verbose:
        for i, (result, _time, command) in enumerate(zip(results, times, test_list)):
            if result == "NG":
                print(f"\r({i+1}/{num_test} {test.name})[{red_str('NG')}] time {_time}s  {command}")
            else:
                print(f"\r({i+1}/{num_test} {test.name})[{green_str('OK')}] time {_time}s  {command}")

#
# utils
#
def red_str(s):
    return "\033[31m" + s + "\033[0m"

def green_str(s):
    return "\033[32m" + s + "\033[0m"

def cyan_str(s):
    return "\033[36m" + s + "\033[0m"


def check_args(args):
    assert args.number_of_parallels >= 4, f"-np or --number_of_parallels must be greater than 3, but got {args.number_of_parallels}"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-o", "--only_commands",
        action="store_true",
        help="display only test commands (no execute test)"
    )
    parser.add_argument(
        "-ns", "--no_test_shared",
        action="store_true",
        help="no test for thread parallel version",
    )
    parser.add_argument(
        "-nd", "--no_test_distributed",
        action="store_true",
        help="no test for mpi parallel version",
    )
    parser.add_argument(
        "-l", "--long",
        action="store_true",
        help="Test by executing 5 minutes (default is 30 seconds)",
    )
    parser.add_argument(
        "-np", "--number_of_parallels",
        type=int,
        default=4,
        help="number of parallels",
    )
    parser.add_argument(
        "-i", "--instance",
        type=str,
        default="./storage/sample_mats/dim80.txt",
        help="instance text file for test",
    )
    parser.add_argument(
        "-v", "--verbose",
        action="store_true",
        help="display standard output for each test",
    )
    args = parser.parse_args()
    check_args(args)
    main(args)
