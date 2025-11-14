import marimo

__generated_with = "0.17.7"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _():
    import marimo as mo
    return (mo,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Tutorial 1: Introduction
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Question 1: Data Types

    Now it's about basic datatypes in python. We will work with the following data
    types:

    * Numeric: integer (int), floating point number (float)

    * Strings: str; e.g. example_string = 'test'

    * Sequence types: list, tuple, range; e.g. example_list = [1,'1',1.0]

    * Mapping data types: dictionary (dict); e.g. example_dict = {1:'test', someList:[1,2,3]}

    * Boolean data type: bool; e.g. example_bool = True

    * Set data type: set


    ---

    Tasks:

    0. Add a comment (#) to every line of code explaining what the code is doing

    1. Go to [UniProt](https://www.uniprot.org/) and retrieve the sequence
    of RNASE1 for Horse and the Common Mink Whale. Copy both and store them
    as strings in two variabels (e.g. Horse_RNASE1 and Whale_RNASE1)

    2. Create a variable "first_alanine" that has an integer value assigned
    referring to the position of the first alanine (A) in the Horse RNASE1. Use it
    to print out the first alanine from the sequence with print(string[first_alanine])

    3. Now hand the same value as a float to the "first_alanine" variable. Try to
    print it again and check the difference.

    4. There is a binding site between position 41 and 45 (KPVNT) in the Whale
    RNASE1 sequence. Print this sniplet from the Horse_RNASE1 variable using the
    [...] notation. Hint: You will need to use a ":" .

    5. Not everyone can understand the amino acid codes. Create a dictionary
    that can be used as a lookup for Lysine, Proline, Valine, Asparagine, and Threonine. It should take the one letter code as a key to return the full name.
    Now print out the AA name of the active site amino acids from task 4 using the
    dictionary.

    6. You want to know how many different amino acids are in the Horse RNASE1.
    Try to count the unique elements by converting the "Horse_RNASE1" into a set in
    combination with the "[len](https://www.w3schools.com/python/ref_func_len.asp)()" function.



    Hint: check https://www.digitalocean.com/community/tutorials/python-data-types
    for help. The web is full of python guides. If you get stuck just google your
    problem. If that still does not help ask the machine overlord ChatGTP. It is
    good practice to find solutions in forums such as stack overflow. ChatGTP will
    not solve all problems (yet).
    """)
    return


@app.cell
def _():
    # 1. Create two variables storing the sequences of RNASE1 horse + whale

    Horse_RNASE1 = 'KESPAMKFERQHMDSGSTSSSNPTYCNQMMKRRNMTQGWCKPVNTFVHEPLADVQAICLQKNITCKNGQSNCYQSSSSMHITDCRLTSGSKYPNCAYQTSQKERHIIVACEGNPYVPVHFDASVEVST'

    Whale_RNASE1 = 'RESPAMKFQRQHMDSGNSPGNNPNYCNQMMMRRKMTQGRCKPVNTFVHESLEDVKAVCSQKNVLCKNGRTNCYESNSTMHITDCRQTGSSKYPNCAYKTSQKEKHIIVACEGNPYVPVHFDNSV'
    return Horse_RNASE1, Whale_RNASE1


@app.cell
def _(Horse_RNASE1):
    # 2. Create an index variable (int)

    first_alanine = Horse_RNASE1.find('A')

    # use it to access string data
    print(f"First alanine on position {first_alanine} -> {Horse_RNASE1[first_alanine]}")
    return


@app.cell
def _(Horse_RNASE1):
    # 3. Create an index variable (float)
    first_alanine_float = float(Horse_RNASE1.find('A'))

    # use it to access string data
    print(f"First alanine on position as float {first_alanine_float} - can not access the str idx")
    return


@app.cell
def _(Whale_RNASE1):
    # 4. Print out active site of the Whale_RNASE1
    print(Whale_RNASE1[40:45])
    return


@app.cell
def _(Whale_RNASE1):
    # 5. Finish the lookup for the aa code

    aa_lookup = {'K':'Lysine', 'P':'Proline', 'V':'Valine', 'N':'Asparagine', 'T':'Threonine'}

    # print the translation
    amino_substring = Whale_RNASE1[40:45]
    print([aa_lookup[protein] for protein in amino_substring])
    return (aa_lookup,)


@app.cell
def _(Horse_RNASE1):
    # 6. Counting unique elements in thr Horse RNASE1
    len(set(Horse_RNASE1))
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Question 2: Loops

    Loops are a basic building block in every programming language and are used
    to automate repetitive tasks.

    The `range()` function is used to generate a sequence of numbers. The syntax of the `range()` function is as follows:

    Syntax: `range([start,] stop [,step])`
    * Start (optional): the starting point of the sequence. It defaults to 0.
    * Stop (required): the endpoint of the sequence. This item will not be included in the sequence.
    * Step (optional): the step size of the sequence. It defaults to 1.

    For example, `range(1,6)` will produce a sequence of numbers from 1 to 5 (be careful, it is not the values of 1 to 6). If one changes it to `range(6)`, the value 0 will also be included because the default start is 0.

    Tasks:

    1. Redo question 5 of the last part. Now use a for loop with the range function
    such that you dont have to copy the same line over and over again.

    2. Print only specific characters from a string. Write a loop that iterates
    the whole Horse RNASE1 sequence but prints only "K"s. Hint: Check out what an
    if statement is and how the == operator can be used for this task.

    3. Adapt the loop from two such that you can count how often a character (and here specifically "K") is
    occuring. Hint: You have to create a new integer variable. You can use the += operator to increase its value within the loop.

    4. Reverse String: Given a string, create a for loop that prints the characters of the string in reverse order. Use the Whale RNASE1 sequence. Hint: Check what happens if you add two strings.
    e.g. "a" + "b"; or "\" + "abc"

    5. (Harder one) Use your knowledge about for loops, if statements, and counting to iterate
    over both sequences at once counting perfect matches. Hint: Check what the "zip" function is doing and how it can be used to create such loop.
    If this is still too boring try to print both Sequences as in the lecture with *
    denoting matches. e.g.:

    ```
    RESP
     ***
    KESP
    ```

    Hint: Check https://www.digitalocean.com/community/tutorials/python-for-loop-example for help. Or Google... or ChatGTP.
    """)
    return


@app.cell
def _(Whale_RNASE1, aa_lookup):
    # 1. Iterate over the positions 40 - 44 while printing from the dictionary in
    # each iteration

    for i in range(40, 45):
        print(aa_lookup[Whale_RNASE1[i]], end=" ")
    return


@app.cell
def _(Whale_RNASE1):
    # 2. Print only K while iterating the whole string

    for char in Whale_RNASE1:
        if char == "K":
            print(char, end=' ')
    return


@app.cell
def _(Whale_RNASE1):
    # 3. Counting K in the Horse RNASE 1
    count = 0
    for char_k in Whale_RNASE1:
        if char_k == "K":
            count += 1
    print(count)
    return


@app.cell
def _(Whale_RNASE1):
    # 4. Reverse string
    reversed_string = ""  # Initialize empty string for reversed characters
    for char_idx in range(len(Whale_RNASE1)-1, 0, -1):
        reversed_string += Whale_RNASE1[char_idx]
    print(reversed_string)
    return


@app.cell
def _(Horse_RNASE1, Whale_RNASE1):
    # 5. Perfect matches
    matches = 0
    for a, b in zip(Whale_RNASE1, Horse_RNASE1):
        matches += a == b
    print(matches)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # Question 3: Functions


    Some basics about functions in python. Functions are defined by the "def"
    statement followed by closed brackets "()" that can contain input variables.
    Afterwards ":" starts the function body which is a block statement and will
    be executed when the function is called. If the function should return something
    the body is closed by the optional return statement followed by whatever
    the output should be. Functions are typically used if complex operation
    have to be applied in multiple parts of the code (there are more reasons to use functions but this is relevant only for advanced programming). Further they are the foundation of object orientated programming.

    Example 1 "No input no return"


    ```
    def show_my_string():
      my_string = "abcd"
      print(my_string)

    show_my_string()

    >>> abcd

    # This is the same as:

    def show_my_string():
      print("abcd")

    ```

    Example 2 "Input and return"



    ```
    def divide(a,b):
      c = a / b

      return c

    result = divide(10, 5)

    print(result)

    >>> 2

    # This is the same as:

    def divide(a,b):
      return a / b

    ```

    ---

    1. Create a function that takes a string variable as input and prints its value. Use it to print the Horse_RNASE1 or Whale_RNASE1

    2. Create a function that counts perfect matches of two strings.

    3. Create a function that takes two strings as input and combines both to a single string. Create two variables one containing the first half of the
    Whale RNASE1 and the other containing the second half of the Horse RNASE1. Use your new function to combine both halfs into a third variable.

    4. How many perfect matches does the newly created RNASE1 have with the Whale RNASE1 and with the Horse RNASE1. Use the function created in task 2.

    5. Now that you have calculated the perfect matches for 3 combinations (Whale vs. Horse, Whale vs. HalfHalf, Horse vs. HalfHalf) you want to compute the average number of perfect matches. Create a function that takes a list of integers or floats and returns the arithmetic mean.

    Hint: Check https://www.digitalocean.com/community/tutorials/how-to-define-functions-in-python-3
    """)
    return


@app.function
# 1. Print a string
def print_string(some_string):
    print(some_string)


@app.cell
def _(Horse_RNASE1):
    print_string(Horse_RNASE1)
    return


@app.function
# 2. Matches
def perfect_matches(string_a, string_b):
    match_count = 0
    for first_val, second_val in zip(string_a, string_b):
        match_count += first_val == second_val
    return match_count


@app.function
# 3. Combine strings
def combine_strings(string_a, string_b):
    return string_a + string_b


@app.cell
def _(Horse_RNASE1, Whale_RNASE1):
    # 4. Find perfect matches
    new_rnase = Whale_RNASE1[:len(Whale_RNASE1)//2] + Horse_RNASE1[len(Horse_RNASE1)//2:]

    print(perfect_matches(new_rnase, Whale_RNASE1))
    print(perfect_matches(new_rnase, Horse_RNASE1))
    return (new_rnase,)


@app.function
# 5. Compute the mean of a list
def comp_mean(numbers):
    sum_val = 0
    for n in numbers:
        sum_val += n
    return sum_val / len(numbers)


@app.cell
def _(Horse_RNASE1, Whale_RNASE1, new_rnase):
    print(comp_mean([perfect_matches(new_rnase, Whale_RNASE1), perfect_matches(new_rnase, Horse_RNASE1), perfect_matches(Whale_RNASE1, Horse_RNASE1)]))
    return


if __name__ == "__main__":
    app.run()
