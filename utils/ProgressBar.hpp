#pragma once
#include <cstdio>
#include <ctime>
#include <functional>
#include <algorithm>
#include <string>
#include <sys/ioctl.h>
#include <unistd.h>

namespace ProgressBar {
  using namespace std;
  function<string()> default_info_func = []() { return string(""); };
  class ProgressBar {
    private:
      int current, total;
      const char *bit = "━";
      const int width = 40;
      time_t start_time;
      char refreshCode;
      function<string()> info_func;
      int console_width, console_height;
    public:
      ProgressBar(int total, char refreshCode='\r', function<string()> info_func=default_info_func) {
        this->total = total;
        this->current = 0;
        this->refreshCode = refreshCode;
        this->info_func = info_func;
        time(&(this->start_time)); 
        struct winsize size;
        ioctl(STDOUT_FILENO, TIOCGWINSZ, &size);
        this->console_width = size.ws_col;
        this->console_height = size.ws_row;
      }
      ~ProgressBar() {}
      void update_info() {
        // 更新一下信息
        this->next(0);
      }
      void next(int update_current = 1) {
        printf("\033[2K"); // 清空当前行
        time_t curtime; time(&curtime);
        double time_used = curtime - (this->start_time);
        if(update_current) {
          current = current + 1; 
        }
        printf("\033[0m\033[1;34m"); // blue
        printf("[% 7.2f%% | % 6.2fs used | % 6.2fh(% 6.2fmin) need | % 4d/%d ]", 
            (double)current/total*100, 
            time_used,
            time_used/current*(total-current)/60/60,
            time_used/current*(total-current)/60,
            current, total);
        printf("\033[0m ");
        int num = (int)((double)current*this->width/total);
        printf("\033[0m\033[1;32m"); // green
        int count_green = max(0, num);
        for(int i = 1, flag = 1 ; i <= this->width ; ++ i) {
          if(count_green > 0) {
            count_green --;
          } else if(flag) {
              flag = 0;
              printf("\033[0m");
              printf(" ");
              printf("\033[0m\033[1;31m"); // red
          }
          printf("%s", bit);
        }
        printf("\033[0m");
        if(current >= total) {
          printf("\033[0m\033[1;32m%s \033[0m\033[1;32m%s\033[0m", this->bit, "Done!     \n");
        } else {
          printf(" \033[0m\033[1;31m%s\033[0m", "Waiting...");
          printf(" [%s] ", this->info_func().c_str());
          putchar(this->refreshCode);
        }
        fflush(stdout);
      }
  };
#define ProgressIteration(iter_var, total_iter, code, info_func) \
  if(int break_flag = 1) \
  for(ProgressBar :: ProgressBar __my_progressbar(total_iter, code, info_func); break_flag ; break_flag=0) \
  for(int iter_var = 0 ; (iter_var) < (total_iter) ; __my_progressbar.next(), ++ (iter_var))
}

