# -*- coding: utf-8 -*-
__author__ = 'Wentong Zhou'

import os, sys, pathlib
from PyQt6 import QtGui, QtCore, QtWidgets, uic
from PyQt6.QtWidgets import (QApplication, QMainWindow, QLabel, QMessageBox,
                             QFileDialog, QGraphicsPixmapItem, QGraphicsScene, QInputDialog, QDialog,
                             QListView, QAbstractItemView, QTreeView, QWidget, QLayout, QVBoxLayout,
                             QHBoxLayout, QGridLayout,QTextEdit, QSpinBox, QAbstractSpinBox,
                             QPushButton, QToolButton, QRadioButton, QCheckBox, QLineEdit, QDoubleSpinBox,
                             QTableWidgetItem, QFrame, QSpacerItem, QSizePolicy, QTableWidget)
from PyQt6.QtGui import (QPixmap, QColor, QPainter, QPen, QFont, QDropEvent, QIcon, QTextCursor,
                         QScreen, QKeyEvent, QTextCharFormat, QSyntaxHighlighter)
from PyQt6.QtCore import QPoint, QTimer, QMimeData, QSize, pyqtSignal, QProcess
from PyQt6.QtCore import Qt as QtCore_Qt


Qt_Keys = QtCore_Qt.Key
Qt_Colors = QtCore_Qt.GlobalColor

QAspectRatioMode = QtCore_Qt.AspectRatioMode
QKeepAspectRatio = QAspectRatioMode.KeepAspectRatio

QAlignmentFlag = QtCore_Qt.AlignmentFlag
QAlignCenter = QAlignmentFlag.AlignCenter

QTransformationMode = QtCore_Qt.TransformationMode
QSmoothTransformation = QTransformationMode.SmoothTransformation

QMessageBox_Abort = QMessageBox.StandardButton.Abort
QMessageBox_Cancel = QMessageBox.StandardButton.Cancel
QMessageBox_Close = QMessageBox.StandardButton.Close
QMessageBox_Discard = QMessageBox.StandardButton.Discard
QMessageBox_Ignore = QMessageBox.StandardButton.Ignore
QMessageBox_No = QMessageBox.StandardButton.No
QMessageBox_NoToAll = QMessageBox.StandardButton.NoToAll
QMessageBox_Ok = QMessageBox.StandardButton.Ok
QMessageBox_Save = QMessageBox.StandardButton.Save
QMessageBox_SaveAll = QMessageBox.StandardButton.SaveAll
QMessageBox_Yes = QMessageBox.StandardButton.Yes
QMessageBox_YesToAll = QMessageBox.StandardButton.YesToAll

QTextCursor_End = QTextCursor.MoveOperation.End


def get_open_directories():
    if not QApplication.instance():
        QApplication(sys.argv)

    file_dialog = QFileDialog()
    file_dialog.setFileMode(QFileDialog.DirectoryOnly)
    file_dialog.setOption(QFileDialog.DontUseNativeDialog, True)
    file_view = file_dialog.findChild(QListView, 'listView')

    # to make it possible to select multiple directories:
    if file_view:
        file_view.setSelectionMode(QAbstractItemView.MultiSelection)
    f_tree_view = file_dialog.findChild(QTreeView)
    if f_tree_view:
        f_tree_view.setSelectionMode(QAbstractItemView.MultiSelection)

    if file_dialog.exec():
        return file_dialog.selectedFiles()

    return []


class Qt_Widget_Common_Functions:
    closing = pyqtSignal()

    def center_the_widget(self, activate_window=True):
        frame_geometry = self.frameGeometry()
        screen_center = QtGui.QGuiApplication.primaryScreen().availableGeometry().center()
        frame_geometry.moveCenter(screen_center)
        self.move(frame_geometry.topLeft())
        if activate_window:
            self.window().activateWindow()

    def closeEvent(self, event: QtGui.QCloseEvent):
        # print("Window {} closed".format(self))
        self.closing.emit()
        if hasattr(super(), "closeEvent"):
            return super().closeEvent(event)

    def open_config_file(self):
        self.config = open_config_file()

    def get_config(self, key, absence_return=""):
        return get_config(self.config, key, absence_return)

    # backward compatible
    def load_config(self, key, absence_return=""):
        return self.get_config(key, absence_return)

    def save_config(self):
        save_config(self.config)


def default_signal_for_connection(signal):
    if isinstance(signal, QPushButton) or isinstance(signal, QToolButton) or isinstance(signal, QRadioButton) or \
            isinstance(signal, QCheckBox):
        signal = signal.clicked
    elif isinstance(signal, QLineEdit):
        signal = signal.textChanged
    elif isinstance(signal, QDoubleSpinBox) or isinstance(signal, QSpinBox):
        signal = signal.valueChanged
    return signal


def disconnect_all(signal, slot):
    signal = default_signal_for_connection(signal)
    marker = False
    while not marker:
        try:
            signal.disconnect(slot)
        except Exception as e:  # TODO: determine what's the specific exception?
            # traceback.print_exc()
            # print(e)
            marker = True


def connect_once(signal, slot):
    signal = default_signal_for_connection(signal)
    disconnect_all(signal, slot)
    signal.connect(slot)



def build_fileDialog_filter(allowed_appendix: list, tags=()):
    """

    :param allowed_appendix: a list of list, each group shows together [[xlsx,log,out],[txt,com,gjf]]
    :param tags: list, tag for each group, default ""
    :return: a compiled filter ready for QFileDialog.getOpenFileNames or other similar functions
            e.g. "Input File (*.gjf *.inp *.com *.sdf *.xyz)\n Output File (*.out *.log *.xlsx *.txt)"
    """

    if not tags:
        tags = [""] * len(allowed_appendix)
    else:
        assert len(tags) == len(allowed_appendix)

    ret = ""
    for count, appendix_group in enumerate(allowed_appendix):
        ret += tags[count].strip()
        ret += "(*."
        ret += ' *.'.join(appendix_group)
        ret += ')'
        if count + 1 != len(allowed_appendix):
            ret += '\n'

    return ret

