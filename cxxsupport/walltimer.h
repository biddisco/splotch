#ifndef PLANCK_WALLTIMER_H
#define PLANCK_WALLTIMER_H

#include <string>
#include <map>

class wallTimerSet
  {
  private:
    struct timer
      {
      double t_acc, t_started;
      timer() { t_acc=t_started=0; }
      };

    std::map<std::string,timer> timers;

  public:
    void start(const std::string &name);
    void stop(const std::string &name);
    void reset(const std::string &name);
    double acc(const std::string &name);

    void report() const;
  };

extern wallTimerSet wallTimer;

#endif
