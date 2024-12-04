
// build % cmake .. -DCMAKE_PREFIX_PATH="/opt/homebrew/Cellar/qt/6.7.2_1/lib/cmake"
#include <QApplication>
#include <QPushButton>

int main(int argc, char *argv[]) {
    QApplication app(argc, argv);
    QPushButton button("Hello, Qt!");
    button.resize(200, 100);
    button.show();
    return app.exec();
}
