# My-bachelor-degree
Численное решение обобённого уравнения Шрёдингера 4 порядка с использованием 
псевдоспектрального мметода Фурье с разделённым шагом (на основе полученного 
ранее аналитического решения в виде светлого солитона).
Реализованы различные варианты программы:
- main-NLS: вместо ОНУШ решается НУШ
- mainTests: решение ОНУШ в виде аналитического солитона
- mainSinNoize / mainWhiteNoize / mainAmpNoize: решение ОНУШ с шумом в начальных данных
(мультипликативный шум / аддитивный шум / вариация амплитуды с сохранением интеграла энергии)
- mainCollision: решение ОНУШ с локализованным шумом в виде гауссова пакета
  
## Требования
- fftw3.3 (https://www.fftw.org/)
- openmp (https://www.openmp.org/)
- g++ (или любой совместимый C++11 компилятор)
- GNU Make
- ОС Linux / macOS (на Windows можно собрать через WSL или MinGW)
- Python 3.8+
- Установленные зависимости (указаны в requirements.txt)
  
## Установка

Установить OpenMP и FFTW3.
Собрать проект с makefile.
Как установить зависимости:
- Для Python: `pip install -r requirements.txt`

## Запуск

Как запустить:
- C++: 'make -f makefile'
Чтобы пересобрать проект, введите 'make -f makefile rebuild'
- Python: `jupyter notebook`

## Пример использования
```bash
#path N T A alpha a b pathRe pathIm
./executeTests 2048 10.0 1.0 1.0 0.5 0.5 logTests/preciseRe.csv logTests/preciseIm.csv
```

