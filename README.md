# Thesis
MA Thesis, Titov S
Strange attractor 'n stuff

Эта программа написана с целью изучения свойств странных аттракторов в дискритезованых пространствах. В данном коде реализован странный аттрактор Лоренца, пространство дискретизовано простой решеткой (кубической).

Алгоритм в общих словах

Алгоритм состоит из двух основных частей
1. Получение наиболее полной траектории аттратора при данной дискретизации
2. Изменение парамтеров дискретизации 

Подробнее о первоой части 

Очевидно, что траектория странного аттрактора перестадет быть бесконечной в дискретизованом пространстве - она распадается на циклы. В нашем случае, мы можем переформулировать задачу из "нахождения наиболее полной траектории" в "нахождение полного набора циклов" странного аттрактора. 

Для поиска отдельного цикла используется Метод Рунге-Кутта 4 порядка для вычисление траектории. Для обнаружения наиболее полного количества циклов, на которые распадается аттрактор мы итеративно вычисляем траектории каждый раз меняя начальные координаты. 

Смена начальных координат для вычисления полного набора проходила по прямой, которая (в случае аттрактора Лоренца) проходила через устойчивые точки аттрактора. Длина прямой превышает выходит за объем аттрактора, что бы гарантировать наибольшее количество наденных циклов. 

Вычисляя каждый цикл, мы сохраняем его длину и количество точек из которых он состоит. 

Подробнее о второй части

В подобном исследовании было обанружено, что при линейном изменении размера решетки дискретизации были выявлены плато циклов на которых наборы оставались идентичными. В данной рабое мы решили сосредоточится на более точном обануржении точек разрыва между плато. Для выполнения этой задачи мы разделили изменение рзмера решетки дискретизации на два типа: 

  Линейный - последовтельное линейное изменение 
  Поиск точки разрыва - когда два сосдених набоа циклов различаются, то мы с помощью метода деления отрезка пополам ищем точку разрыва до   опредленной точности

Пользуясь двумя этими методами мы можем получить довольно точный потретрет поведения наборов циклов 

Теперь о программе

В програме реализовано три класса 

1. Lorenz - класс содержащий все переменные и фукнции касающиеся системы лоренца
2. Grid - класс решетки - ее размер, шаг и прочее 
3. Stat - отвчает за сбор статистики с полученных раекторий (на данныыый момент подсчет длинны цикла)

В рамках этих классов происходит вся работа программы - с более подробным описанием можно ознакомится в комментария в коде

Детали

Хранение и вычисление траеторий происходит в одномерном массиве, где tr[0], tr[1], tr[2] соответсвенно x y и z координаты точки аттрактора. Именно из за этого в коде можно встретить иногда множители при обращении к массиву

Поиск цикла реализован с помощью хэш-массива. Три координаты каждой точки возвращают уникальное значение хэш функции. При одинаковых значениях шэшфункций производится сравнение координат. Если они равны - цикл найден.

Динамическая память. Из за упрощеной модели хранения траекторий динамическое расширение памяти пришлость реалиховывать руками

Степень готовности.

На данный момент программа работает, но не проверялась на массивном вычислении. Из того, что надо доделать - вывод данных в соответсвующем формате
