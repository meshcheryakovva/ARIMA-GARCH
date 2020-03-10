function arimapredict
  % Функция прогнозирования нестационарных временных рядов
  % на основе ARIMA/GARCH моделей.
  
  % Версия для запуска под MATLAB.
  
  % Используются функции пакетов расширения:
  % Econometrics, Statistics, Optimization Toolboxes
  
  % На входе 3 файла:
  % predictparam.txt - файл параметров прогноза
  % learn.txt - исходные данные для построения модели ряда
  % verif.txt - исходные данные для верификации (сравнения с прогнозом)
  
  % на выходе 2 файла:
  % ypredict.txt - файл, содержащий точечные прогнозы (т.е. поведение тикеров в среднем)
  % proposit.txt - файл, содержащий вероятности попадания в доверительный интервал
  clear; clc % очистка командного окна, рабочей области 
  close all  % закрытие графиков
  
  % считываем файл параметров прогноза:
  % 1) количество последних наблюдений, используемых для прогнозирования;
  % 2) количество часов вперед;
  % 3-6) ширина доверительного интервала для Open, High, Low, Close
  load predictparam.txt
  % используем last последних уровней ряда, т.е. last последних часов
  % (если уровни временного ряда пропущены, т.е. нет торгов,
  % пропуски заполняются предыдущими значениями)
  last = predictparam(1); 
  predictions = predictparam(2); % количество часов вперед (прогноз);
  UO = predictparam(3); % ширина доверительного интервала для O
  UH = predictparam(4); % ширина доверительного интервала для H
  UL = predictparam(5); % ширина доверительного интервала для L
  UC = predictparam(6); % ширина доверительного интервала для C
  
  % считываем временные ряды:
  % имя файла, содержащего тикеры
  fileName = 'learn.txt'; % данные для построения модели ряда EUR/USD
  % считываем и заполняем пропущенные строки по времени.
  % на выходе - моменты времени в часах и тикеры
  [Time, O, H, L, C] = readDayOHLC(fileName, last);
    
  % тикеры на момент окончания суток отмечаем отдельно (для наглядности на графике)
  shift = mod(Time(1),24); % сдвиг ряда относительно начала суток
  TimeDays = Time(25-shift:24:end); % моменты окончания суток
  ODays = O(25-shift:24:end); % соответствующие уровни ряда Open
  HDays = H(25-shift:24:end); % соответствующие уровни ряда High
  LDays = L(25-shift:24:end); % соответствующие уровни ряда Low
  CDays = C(25-shift:24:end); % соответствующие уровни ряда Close
  
  % вызов функций формирования моделей рядов,
  % расчет точечных прогнозов и вероятностей попадания прогнозов
  % в доверительный интервал,
  % а также ширины дов. интервала для 99% вероятности попадания
  [yOpen, pOpen, o99]   = ARIMA(Time, O, UO, predictions, TimeDays, ODays, 'Open');
  [yHigh, pHigh, h99]   = ARIMA(Time, H, UH, predictions, TimeDays, HDays, 'High');
  [yLow, pLow, l99]     = ARIMA(Time, L, UL, predictions, TimeDays, LDays, 'Low');
  [yClose, pClose, c99] = ARIMA(Time, C, UC, predictions, TimeDays, CDays, 'Close');
  
  % точечный прогноз
  ypredict = [yOpen yHigh yLow yClose];
  % вероятности попадания в заданный пользователем интервал
  proposit = [pOpen pHigh pLow pClose];
  
  % границы интервального прогноза 99%
  lb = [yOpen-o99 yHigh-h99 yLow-l99 yClose-c99];
  ub = [yOpen+o99 yHigh+h99 yLow+l99 yClose+c99];
  
  % график прогноза
  figure
  plot(1:predictions, ypredict) % точечный прогноз
  grid on
  title('Точечные прогнозы и довер. интервалы (99%). Точки - фактич.наблюдения')
  legend('Open', 'High', 'Low', 'Close')
  hold on
  plot(1:predictions, lb, ':') % интервальный прогноз
  plot(1:predictions, ub, ':')
  
  fileName = 'verif.txt'; % данные для верификации прогноза
  [~, Ov, Hv, Lv, Cv] = readDayOHLC(fileName, predictions);
  plot(1:predictions, Ov, 'b.') % фактические наблюдения
  plot(1:predictions, Hv, 'g.')
  plot(1:predictions, Lv, 'r.')
  plot(1:predictions, Cv, 'm.')
  xlabel('Прогнозные моменты времени (час)')   
  
  % сохраняем результат в текстовые файлы
  % в столбцах - значения для Open, High, Low, Close
  % строки соответствуют прогнозным моментам времени в часах
  % ypredict.txt - файл, содержащий точечные прогнозы
  save ypredict.txt ypredict -ascii 
  % proposit.txt - файл, содержащий вероятности попадания в доверительный интервал
  save proposit.txt proposit -ascii 
end


function [yprediction, propab, U99] = ARIMA(Time, y, U, predictions, TimeDays, yDays, yName)
  % выбор модели, моделирование и прогнозирование ряда y
 
  % выход:
  % yprediction - точечные прогнозы
  % propab - вероятности попадания прогноза в пользовательский интервал
  % U99 - ширина дов. интервала при 99% вероятности прогнозирования (только для верификации)
  
  % вход:
  % Time - моменты времени наблюдений (часы)
  % y - уровни врем. ряда (по часам)
  % U - заданная пользователем ширина доверительного интервала
  % predictions - заданное пользователем количество часов вперед для прогноза
  % TimeDays - моменты времени наблюдений для окончания каждых суток
  % yDays - соответствующие уровния ряда
  % yName - текстовая строка, имя показателя
  
  disp(yName) % вывод в командное окно
  n = length(y); % длина временного ряда
  % Графики рядов
  figure
  plot(Time, y, '.b')
  hold on
  plot(TimeDays, yDays, '.r')
  title(yName); xlabel('Время, ч')
  grid on
  
  % Строим модель нестационарного временного ряда: ARIMA(p, d, q)
  
  % переход к стационарному ряду путем взятия разностей
  d = 1; % степень интегрируемости
  % Расширенный тест Дики-Фуллера, уровень ошибки 0.1
  % модель авторегрессии со смещением (тренда нет, ненулевое матожидание)
  h = adftest(y, 'model', 'ARD', 'alpha', 0.1);
  
  while ~h, % проверка нуль-гипотезы
      dy = diff(y, d); % разностный ряд степени d
      % Расширенный тест Дики-Фуллера, уровень ошибки 0.1
      % модель авторегрессии (тренда нет, нулевое матожидание разностного ряда)
      h = adftest(dy, 'model', 'AR', 'alpha', 0.1);
      d = d + 1; % увеличение степени интегрируемости
  end
  dy = diff(y, d); % разностный стационарный ряд степени d
  fprintf('Порядок разностного ряда d = %d\n', d)
  
  % начальные условия для интегрирования прогнозного ряда
  iy = zeros(d,1);
  iy(1) = y(1); % запоминаем первый уровень исх. ряда
  for i = 1:d-1, % запоминаем первый уровень разностного ряда каждой степени
      dify = diff(y, i);
      iy(i+1) = dify(1);
  end
  
  % строим ARMA(p, q) для стационарного разностного ряда dy
  % подбор порядков модели
  Semin = std(y); % ищем минимальное стандартное отклонение остатков
  for p = 1:4,     % параметр AR - размерность авторегрессионной части модели
      for q = 0:3; % параметр МА - размерность скользящ.среднего в составе модели
          % задаем структуру модели (спецификация):
          % предполагаем "толстые хвосты" функции плотности вероятности,
          % выбираем t-распределение
          SpecDy = garchset('R', p, 'M', q, 'Display', 'Off', 'Distribution', 'T');
          % подгонка ARMA:
          % EstSpecDy - структура, содержащая параметры модели
          % Innov - инновации (остатки)
          [EstSpecDy, ~, ~, Innov] = garchfit(SpecDy, dy);
          
          % расчетные значения ARIMA находим интегрированием:
          yhat = dy - Innov; % расчетные значения (ARMA) разностного ряда
          for i = d-1 : -1 : 0, % интегрирование в цикле
              dyhat = yhat; % это оценки приращений
              if i==0,
                  dyi = y;
              else
                  dyi = diff(y, i); % предыдущий ряд разностей
              end
              % оценка разностей
              dyhat = dyi(1:end-1) + dyhat; % прибавляем оценки приращений
              % оценка исходного ряда
              yhat = [iy(i+1); dyhat]; % добавляем начальные условия
          end
          
          % ряд остатков
          e = y - yhat; % остатки модели ARIMA
          Se = std(e) ; % СКО
          % ищем параметры ARMA, соответствующие минимальному СКО ошибки
          if Se < Semin,
              pbest = p; % лучшая размерность модели
              qbest = q;
              Semin = Se; % минимальное стандартное отклонение остатков
              ESTSPECDY = EstSpecDy; % лучшие параметры модели ARMA
              YHAT = yhat; % расчетные значения лучшей модели ARIMA
              E = e; % лучшие остатки ARIMA
          end
      end
  end
  
  yhat = YHAT; % расчетные значения ARIMA для лучшей модели ряда y
  e = E; % соответствующие остатки
  % точечный прогноз (ARMA) разностного ряда
  [~, meanForecast] = garchpred(ESTSPECDY, dy, predictions);
  
  fprintf('Наилучшая модель - ARIMA(%d, %d, %d)\n', pbest, d, qbest)
  fprintf('Средняя относительная ошибка %f\n', mean(abs(e)./y))
  
  % В операторном виде находим обобщенный нестационарный оператор
  % авторегрессии ARIMA, затем находим функцию реакции на импульсы
  % (для прогнозирования дисперсии остатков на несколько шагов вперед)
  ARcell = num2cell([1 -ESTSPECDY.AR]); % коэффициенты авторегрессии в структуре ARMA
  AL = LagOp(ARcell); % операторный полином (левая часть ARMA)
  MAcell = num2cell([1 ESTSPECDY.MA]); % коэффициенты скольз.среднего в структуре ARMA
  MA = LagOp(MAcell); % операторный полином (правая часть ARMA)
  I1 = LagOp({1 -1}); % разностный оператор
  for i = 1:d, % ARMA приводим к ARIMA
      AL = AL * I1;
  end
  warning off % предупреждения обусловлены нестационарностью ряда
  % приводим ARIMA к виду MA:
  % полиномиальное деление в операторном виде
  PSIop = mldivide(AL, MA, 'Degree', predictions-1); 
  warning on
  % извлекаем коэффициенты
  PSI = cell2mat(toCellArray(PSIop)); % функция реакции на импульсы
  % множитель, учитыв. рост СКО прогнозирования
  % при удалении от известных данных:
  SeGrow = cumsum(PSI .^ 2) .^ 0.5;   
  
  % точечный прогноз ARIMA:
  % интегрируем d раз разностный ряд с прогнозом ARMA
  pred = cumsum([iy(d); dy; meanForecast]);
  for i = 1:d-1,
      pred = cumsum([iy(d-i); pred]);
  end
  % pred(n+1:end) - это точечный прогноз по ARIMA без учета гетероскедастичности
  
  % Оценка волатильности с помощью дисперсии остатков.
  % строим GARCH-модель для ряда остатков
  p = pbest;
  q = qbest;
  % тест Ингла на наличие гетероскедастичности остатков, уровень ошибки 0.1
  h = archtest(e, 'lags', p + q, 'alpha', 0.1);
  if ~h, % гетероскедастичность не обнаружена
      yprediction = pred(n+1:end); % точечный прогноз
      MeanTime = mean(Time); % среднее значение времени наблюдений
      % стандартная ошибка прогнозирования c учетом дальности прогноза
      Sigma = sqrt(sum((e - mean(e)).^2)/(n-2)) * (1 + 1./n + ...
          ((((n+1):(n+predictions))' - MeanTime).^2 ./ sum((Time - MeanTime).^2))).^0.5;
      Sigma = SeGrow' .* Sigma; % учет нестационарности ARIMA
      disp('гетероскедастичность остатков не обнаружена')
  else % предполагается гетероскедастичность
      % структура модели GARCH
      Spec = garchset('R', p+q, 'M', 0, 'VarianceModel', 'GARCH', ...
                      'Display', 'Off', 'Distribution', 'T');
      % оценка параметров модели
      EstSpec = garchfit(Spec, e);
      % прогнозирование СКО и среднего для ряда остатков
      [~, meanForecast, ~, MeanRMSE] = garchpred(EstSpec, e, predictions);
      % учет гетероскедастичности при точечном прогнозе - коррекция
      yprediction = pred(n+1:end) + meanForecast; % скорректированный точечный прогноз
      disp('гетероскедастичность остатков обнаружена')
      Sigma = SeGrow' .* MeanRMSE; % условное среднеквадратическое отклонение
  end
  t = U ./ Sigma; % квантили, соотв. заданной пользователем ширине дов. интервала
  propab = 2 * tcdf(t, n-2) - 1; % вероятность попадания в доверительный интервал
  
  p99 = 0.01; % задаем вероятность ошибки для верификации
  U99 = tinv(1 - p99/2, n-2) .* Sigma; % ширина дов. интервала для вероятности 99%
 
  % график прогноза ARIMA:
  prdc = int32(predictions);
  AllTime = Time(1)+(0:n-1 + prdc);
  plot(AllTime(n+1:end), pred(n+1:end), '.g-')
  % ряд расчетных и прогнозных значений
  yhatpred = [yhat; yprediction];
  plot(AllTime, yhatpred, '-m.')
  % график доверительного интервала пользователя
  errorbar(AllTime(n+1:end), yprediction, repmat(U, predictions, 1), 'sk')
  % подписи - вероятности попадания в интервал
  for i=1:predictions,
      text(double(AllTime(n+i)), yprediction(i)+2*U, num2str(propab(i)), 'FontSize', 10)
  end
  legend('y (час)', 'y (сутки)', 'прогноз ARIMA', 'y расчетн. и точ. прогноз GARCH', 'дов. интервал')
  
  % график остатков
  figure;
  plot(Time, e, '.b-')
  grid on
  title('Остатки - разница между y и ARIMA')
  xlabel('Время, ч')
  
  % прогнозирование разностного ряда - график
  figure
  plot(dy, '.b-')
  hold on
  plot(dy-Innov,'.r-')
  grid on
  title(strcat('моделирование разностного ряда \Delta^', num2str(d), 'y'))
  xlabel('Время, ч')
  legend('dy', 'ARMA(p, q)')
end

function [T, rO, rH, rL, rC] = readDayOHLC(HFileName, oflast)
  % Считывание тикеров из файла часовых цен
  % HFileName - строка, содержащая путь к файлу
  % oflast - количество используемых последних уровней ряда
  
  fid = fopen(HFileName);
  % читаем заголовки столбцов
  textscan(fid, '%s', 7, 'delimiter', '>');
  % читаем строки данных
  Cstr = textscan(fid, '%*s %d %d %f %f %f %f', 'delimiter', ',');
  fclose(fid);
  
  % количество наблюдений - длина ряда
  nmin = int32(min(length(Cstr{:,1}), oflast) - 1);
  rDates = Cstr{:, 1}(end-nmin:end); % даты
  rHours = round(Cstr{:, 2}(end-nmin:end) ./ 10000); % часы (от 01 до 24)
  
  CD = [Cstr{:, 3:6}]; % выделение массива чисел из массива ячеек
  CDred = CD(end-nmin:end, 1:4); % берем последние last строк
  
  % вызов функции заполнения пропущенных наблюдений
  [T, OHLC] = ReplicateMissed(rDates, rHours, CDred);
  T  = T(end-nmin:end); % время (в часах без пропусков)
  rO = OHLC(end-nmin:end, 1); % часовые цены Open
  rH = OHLC(end-nmin:end, 2); % часовые цены High
  rL = OHLC(end-nmin:end, 3); % часовые цены Low
  rC = OHLC(end-nmin:end, 4); % часовые цены Close
  
end

function [t, augtseries] = ReplicateMissed(days, hours, tseries)
  % функция повтора пропущенных значений временного ряда
  % если уровень временного ряда отсутствует, вставляется предыдущее значение

  % t - время с равномерным шагом
  % augtseries - дополненный многомерный временной ряд без пропущенных уровней

  % days  - массив дат
  % hours - массив часов
  % tseries - многомерный временной ряд с пропущенными наблюдениями
  
  ddays = diff(days); % разница дат
  changd = find(ddays ~= 0) + 1; % ищем моменты смены дат
  % приводим часы к монотонной последовательности
  for j = 1:length(changd)
    hours(changd(j):end) = hours(changd(j):end) + 24;
  end
  dhours = diff(hours) - 1; % количество пропущенных часов
  
  missed = find(dhours > 0); 
  n = length(missed);
  htimes = missed;
  dhours = double(dhours);
  for j = 1:n,
      position = missed(j); % куда копируем
      whattocopy = tseries(position,:); % что копируем
      howmtimes = dhours(htimes(j)); % сколько раз
      tseries = [tseries(1:position, :); ...
                 repmat(whattocopy, howmtimes, 1); ...
                 tseries(position+1:end, :)];
      missed = missed + howmtimes;
  end
  t = hours(1):hours(end);
  augtseries = tseries;
end
