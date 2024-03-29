#pragma once
#ifndef UTILS_H
#define UTILS_H

#include <chrono>
#include <fstream>
#include <iterator>
#include <mutex>
#include <random>
#include <regex>
#include <string>
#include <sstream>
#include <thread>
#include <vector>

//===============================================================//
// System
//===============================================================//
void debug(bool condition, std::string info, int status = -1) {
    if (condition) {
        if (status == -1) {
            std::cerr << "warning: " << info << "\n";
        } else {
            std::cerr << "error: " << info << "\n";
            exit(status);
        }
    }
}

//===============================================================//
// Random
//===============================================================//
std::random_device rd;
std::mt19937 e2(rd());
std::uniform_real_distribution<> E(0,1);

//===============================================================//
// Threading
//===============================================================//
std::mutex m;

class threader {
private:
    int N_THREADS;
public:
    std::vector<std::thread> THREADS;

    threader() {
        N_THREADS = std::thread::hardware_concurrency() - 2;
        THREADS   = std::vector<std::thread>(N_THREADS);
    }

    threader(int N_THREADS) : N_THREADS(N_THREADS) {
        THREADS   = std::vector<std::thread>(N_THREADS);
    }

    void join() {
        for (auto& th : THREADS) {
            th.join();
        }
    }

    void lock() {
        m.lock();
    }

    void unlock() {
        m.unlock();
    }
};

//===============================================================//
// Number utils
//===============================================================//
template <typename T>
int sign(T var) {
    if (var < 0) return -1;
    else return 1;
}

//===============================================================//
// String utils
//===============================================================//

/**
 * @brief Tokenize a string with a given delimiter. It uses a blank space by default.
 * 
 * @param s String to tokenize.
 * @param delimiter 
 * @return std::vector<std::string> 
 */
std::vector<std::string> tokenize(std::string s, std::string delimiter = " ") {
    size_t pos = 0;
    std::vector<std::string> tokens;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        tokens.push_back(s.substr(0, pos));
        s.erase(0, pos + delimiter.length());
    }
    tokens.push_back(s);
    return tokens;
}

/**
 * @brief Tabulates the string.
 * 
 * @param tab Tabulation type.
 * @param n Indentation.
 * @return std::string 
 */
std::string bleeding(std::string tab, int n) {
    std::string tabulated = "";
    for (int i = 0; i < n; i++) {
        tabulated += tab;
    }
    return tabulated;
}

//===============================================================//
// Vector utils
//===============================================================//

// Print vector.
template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<T>& v) {
    for (auto& i : v) {
        os << i << " ";
    }
    return os;
}

// Print matrix.
template <typename T>
std::ostream& operator<<(std::ostream& os, std::vector<std::vector<T>>& v) {
    for (auto& i : v) {
        for (auto& j : i) {
            os << j << " ";
        }
        os << "\n";
    }
    return os;
}

// Insert element in sorted vector and return its correspondent index.
template <typename T>
int insert(std::vector<T> &v, T value ) {
    auto at = v.size();
    for (at = 0; at < v.size() && (v[at] <= value); at++);
    // typename std::vector<T>::iterator it = std::lower_bound(
    //     v.begin(), v.end(), value, std::less_equal<T>()
    // ); // Find proper position in descending order.
    v.insert(v.begin() + at, value);
    //v.insert(it, value); // Insert before iterator it.
    return at;
}


//===============================================================//
// Regex utils
//===============================================================//

const std::string int_d = "[^0-9]";    // All digits.
const std::string float_d = "[^0-9.]"; // All digits and points.
const std::string new_line = "\\r";    // Special scape character \r.

// Getline + regex, obtiene la linea y la limpia.
std::string get_line(std::ifstream& in, const std::string r = "\\r") {
    std::string s;
    getline(in, s);
    return std::regex_replace(s, std::regex(r), "");
}

std::string replace(std::string s, const std::string& regex, std::string r = "") {
    return std::regex_replace(s, std::regex(regex), r);
}

//===============================================================//
// File directories utils
//===============================================================//

std::string get_path(std::string path) {
    return path.substr(0, path.find_last_of("/\\") + 1);
}

std::string get_filename(std::string path) {
    return path.substr(path.find_last_of("/\\") + 1);
}

//===============================================================//
// Progress bar
//===============================================================//

// Elapsed time.
std::string elapsed_time(std::chrono::seconds s) {

    auto h = std::chrono::duration_cast<std::chrono::hours>  (s); s -= h;
    auto m = std::chrono::duration_cast<std::chrono::minutes>(s); s -= m;

    std::string result("");
    if (h < std::chrono::hours(10))   result += "0";
    result += std::to_string(h/std::chrono::hours(1)) + ":";
    if (m < std::chrono::minutes(10)) result += "0";
    result += std::to_string(m/std::chrono::minutes(1)) + ":";
    if (s < std::chrono::seconds(10)) result += "0";
    return result + std::to_string(s/std::chrono::seconds(1));

}

void flush_stream(std::ostream& os) {
    os.flush();
    os << std::endl;
}

std::chrono::steady_clock::time_point now() {
    return std::chrono::steady_clock::now();
}

using Color = std::string;
const Color orange = "\033[38;2;209;161;22;1m";
const Color white  = "\033[0m";

using Style = std::vector<std::string>;
const Style STYLE1 = {"=",">"," "};
const Style STYLE2 = {"#",">","_"};

class Progress_bar {
private:
    int bar_width;   // Bar width.
    int c_progress;  // Current progress.
    int t_progress;  // Top progress.
                     // Symbols to be drawn.
    bool initialized = false;
    std::vector<std::string> style;
                     // - [0]: solid bar filling.
                     // - [1]: Current position pointer symbol.
                     // - [2]: empty bar filling.
    std::string eta;
                     // Interval time between refresh.
    std::chrono::milliseconds update_t;
    std::chrono::steady_clock::time_point task_start;
    std::chrono::steady_clock::time_point timer;

public:

    Progress_bar() {};
    Progress_bar(int bar_width, std::vector<std::string> style, int update_t) {
        this->bar_width  = bar_width;
        this->style = style;
        this->update_t = std::chrono::milliseconds(update_t);
    }

    void init(std::ostream& os, std::string task, int t_progress) {
        initialized = true;
        os << "[ STARTING " << task << "... ]\n";
        this->c_progress = 0;
        this->t_progress = t_progress;
        this->eta = "NA:NA:NA";
        task_start = timer = now();
    }

    void update(std::ostream& os, std::chrono::nanoseconds time) {

        if (!initialized) {
            std::cerr << "warning: the progress bar has not been initialized.\n";
            return;
        }

        double progress  = ((double) ++c_progress / (double)t_progress);

        auto n = now();
        if (std::chrono::duration_cast<std::chrono::milliseconds>(n - timer) > update_t) {
            timer = n;
            time *= (t_progress - c_progress);
            eta = elapsed_time(std::chrono::duration_cast<std::chrono::seconds>(time));
        }

        int pos = bar_width * progress;
        std::string s_bar("ETA: " + eta + " [");
        for (int i = 0; i < bar_width; ++i) {
            if (i < pos) s_bar += style[0];
            else if (i == pos) s_bar += style[1];
            else s_bar += style[2];
        }
        s_bar += "] " + 
            static_cast<const std::stringstream&> (std::stringstream() << 
                std::fixed << std::setprecision(2) << progress * 100.).str() 
            + "%\r";
        os << s_bar;
        os.flush();
    }

    void kill(std::ostream& os) {
        flush_stream(os);
        os << "FINISHED IN "
            << elapsed_time(std::chrono::duration_cast<std::chrono::seconds>(now() - task_start))
            << std::endl;
        initialized = false;
    }
    
};

#endif