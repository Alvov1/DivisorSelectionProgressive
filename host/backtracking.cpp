#include "HashSelectionHost.h"

namespace HashSelection {
    std::optional<Word> backtracking(const Word& forWord, const Closure& onClosure) {
        const auto& [initPattern, initSize] = forWord;

        MyStack<thrust::tuple<Char, uint8_t, uint8_t>> extensionsStack {};
        MyStack<thrust::pair<Char, int8_t>> permutationsStack {};
        unsigned wordPosition = 0;

        /* Skipping first N non-vowel characters inside pattern. */
        for (; !isVowel(initPattern[wordPosition]) && wordPosition < initSize; ++wordPosition)
            extensionsStack.push({initPattern[wordPosition], 1, 1});

        do {
            if (wordPosition < initSize) {
                /* Count the number of repetition vowels. */
                uint8_t vowelsCount = 1;
                for (unsigned i = wordPosition + 1; isVowel(initPattern[i]) && initPattern[i] == initPattern[wordPosition]; ++vowelsCount, ++i);

                /* Pushing new value in stack */
                extensionsStack.push({
                    initPattern[wordPosition],
                    vowelsCount,
                    (isVowel(initPattern[wordPosition]) && vowelsCount == 1) ? uint8_t(2) : vowelsCount
                                     });
                wordPosition += vowelsCount;

            } else {
                /* Found new pattern. Pushing into buffer. */
                const Word extension = [&extensionsStack] {
                    Word result {}; auto& [data, size] = result;
                    for(unsigned i = 0; i < extensionsStack.size(); ++i)
                        for(unsigned j = 0; j < thrust::get<2>(extensionsStack[i]); ++j)
                            data[size++] = thrust::get<0>(extensionsStack[i]);
                    return result;
                } ();
                const auto value = [&onClosure] (const Word& word, MyStack<thrust::pair<Char, int8_t>>& stack) {
                    const auto& [pattern, patternSize] = word;

                    stack.push({pattern[0], -1});

                    while(!stack.empty()) {
                        if(stack.size() >= patternSize) {
                            const Word permutation = [] (const MyStack<thrust::pair<Char, int8_t>>& stack) {
                                Word word {}; auto& [data, size] = word;
                                for(unsigned i = 0; i < stack.size(); ++i)
                                    data[size++] = stack[i].first;
                                return word;
                            } (stack);
                            if(onClosure(permutation)) return std::optional {permutation };

                            thrust::pair<Char, int8_t> current {}; auto& [sym, varPos] = current;
                            do {
                                current = stack.pop(); ++varPos;
                                const auto& variants = getVariants(pattern[stack.size()]);
                                if(varPos < variants.size) break;
                            } while (!stack.empty());

                            const auto& variants = getVariants(pattern[stack.size()]);
                            if(varPos < variants.size || !stack.empty())
                                stack.push({variants[varPos], varPos});
                        } else
                            stack.push({pattern[stack.size()], -1});
                    }

                    return std::optional<Word> {};
                } (extension, permutationsStack);
                if(value.has_value()) return value;

                /* Popping values from the stack until it's empty or another vowel is found. */
                thrust::tuple<Char, uint8_t, uint8_t> current {};
                do {
                    current = extensionsStack.pop();
                    wordPosition -= thrust::get<1>(current);
                } while (!extensionsStack.empty() && thrust::get<2>(current) < 2);

                if (thrust::get<2>(current)-- > 1)
                    extensionsStack.push(current);
                wordPosition += thrust::get<1>(current);
            }
        } while (!extensionsStack.empty());

        return std::optional<Word> {};
    }

    std::optional<Word> process(const std::vector<Word>& words, const Closure& onClosure) {
        for(unsigned i = 0; i < words.size(); ++i) {
            const auto value = backtracking(words[i], onClosure);
            if(value.has_value()) return value;
//                else Time::cout << "Word " << word << " completed." << Time::endl;

            const auto percent = static_cast<unsigned>(static_cast<double>(i) / static_cast<double>(words.size()) * 100);
            if(percent > 0 && percent % 5 == 0)
                Time::cout << "Completed " << percent << "%." << Time::endl;
        }

        return std::optional<Word> {};
    }
}


std::optional<Word> foundPermutations(const Word& forWord, const std::function<bool(const Word&)>& onClosure){
    const auto& [pattern, patternSize] = forWord;

    std::vector<std::pair<char, short>> stack; stack.reserve(patternSize);
    stack.emplace_back(pattern[0], -1); // Symbol; Amount

    while(!stack.empty()) {
        if(stack.size() >= patternSize) {
            const Word united = [](const std::vector<std::pair<char, short>>& stack) {
                Word word {}; auto& [data, size] = word;
                for(const auto& [sym, _]: stack)
                    data[size++] = sym;
                return word;
            } (stack);
            if(onClosure(united)) return { united };

            unsigned nextPosition = 0;
            do {
                nextPosition = stack.back().second + 1;
                stack.pop_back();

                const auto& variants = getVariants(pattern[stack.size()]);
                if(nextPosition < variants.size()) break;
            } while (!stack.empty());

            const auto& variants = getVariants(pattern[stack.size()]);
            if(nextPosition < variants.size() || !stack.empty())
                stack.emplace_back(variants[nextPosition], nextPosition);
        } else
            stack.emplace_back(pattern[stack.size()], -1);
    }

    return {};
}

bool isReachable(int maze[N][M])
{
    // Initially starting at (0, 0).
    int i = 0, j = 0;

    stack<node> s;

    node temp(i, j);

    s.push(temp);

    while (!s.empty()) {

        temp = s.top();
        int d = temp.dir;
        i = temp.x, j = temp.y;

        temp.dir++;
        s.pop();
        s.push(temp);

        if (i == fx and j == fy) {
            return true;
        }

        if (d == 0) {
            if (i - 1 >= 0 and maze[i - 1][j] and
                visited[i - 1][j]) {
                node temp1(i - 1, j);
                visited[i - 1][j] = false;
                s.push(temp1);
            }
        }

        else if (d == 1) {
            if (j - 1 >= 0 and maze[i][j - 1] and
                visited[i][j - 1]) {
                node temp1(i, j - 1);
                visited[i][j - 1] = false;
                s.push(temp1);
            }
        }

        else if (d == 2) {
            if (i + 1 < n and maze[i + 1][j] and
                visited[i + 1][j]) {
                node temp1(i + 1, j);
                visited[i + 1][j] = false;
                s.push(temp1);
            }
        }

        else if (d == 3) {
            if (j + 1 < m and maze[i][j + 1] and
                visited[i][j + 1]) {
                node temp1(i, j + 1);
                visited[i][j + 1] = false;
                s.push(temp1);
            }
        }

        else {
            visited[temp.x][temp.y] = true;
            s.pop();
        }
    }

    return false;
}